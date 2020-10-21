### Helper classes and functions to create a system
import sim
from collections import namedtuple

from dotdict import dotdict
import forcefield as ff

### wishlist
#1) bonding


### DataTypes (for holding parameters), defined to enable dot notation)
#   python 2.7 doesn't support defaults option for namedtuple...
#   alternative is to use dictionary and key notation
#Topology = namedtuple('Topology','atom_names atom_charges molecule_names molecule_defs bonding',
#        defaults = [ [],{}, [],{}, {} ] )
#        i.e. expected elements are: list, dict, list, dict, dict
Topology = namedtuple('Topology','atom_names atom_charges molecule_names molecule_defs bonding')

def make_sys_params(*args,**kwargs):
    sys_params = dotdict({
        'sys_name':None,
        'num_mols':{},
        'dt':0.05,
        'temp':1.,
        'pressure':0,
        'barostat_options':{'tension':None,'pressure_axis':None},
        'langevin_gamma':1.,
        'steps_equil':None,
        'steps_prod':None,
        'steps_stride':None,
        'md_engine':'openmm',
        'units':None,
        'box':None,
        'traj':None,
        })
    sys_params.update(*args,**kwargs)
    print('system parameters: {}'.format(sys_params))
    return sys_params

def make_special_options(*args,**kwargs):
    special_options = dotdict({
        'sys_name':None,
        'neutralize':False,
        'recalc':False
        })
    special_options.update(*args,**kwargs)
    return special_options

### Convenience functions for defining bonding graphs
def bondgraph_simple(bond_sequence):
    """return simple bonding graph, only for linear/simply branched structures
    
    Notes
    -----
    if bond_sequence is an integer, assume is that many beads linearly bonded
    otherwise, bond_sequence should be a list or tuple, that is the sequence of # beads in each (linearly-bonded) branched unit, i.e. [5,1,5,1,5] for (5-bead) grafts spaced by one monomer
    """
    if isinstance(bond_sequence,int):
        #then this is a simple linear bonding graph
        bond_sequence = [1]*bond_sequence
    bond_list = []    
    previous_backbone_index = 0
    backbone_index = 0
    for num in bond_sequence:
        if backbone_index > 0: #connect to previous backbone monomer/index
            #print("...Bonding {} to {}".format(previous_backbone_index,backbone_index))
            bond_list.append( (previous_backbone_index,backbone_index) )
        for i in range(backbone_index,backbone_index+num-1):
            #print("...Bonding {} to {}".format(i,i+1))
            bond_list.append( (i,i+1) )

        previous_backbone_index = backbone_index
        backbone_index = backbone_index + num

    return bond_list


### Create the system(s)
def create_system(sys_params, topology, interactions, special_options=None, load=True):
    print('=== Creating System ===')
    print('... Initializing Topology ...')
    sim_atom_types = {} #Sys.World.AtomTypes should be same thing, except only keeping atom types that make it into the world, and in list form
    sim_mol_types = [] #Sys.World should be same thing, except only keep mol types that make it into the world, and in list form

    # 2a) define atoms
    for atom_name in topology.atom_names:
        charge = 0. if (special_options and special_options['neutralize']) else topology.atom_charges[atom_name]
        atom_type = sim.chem.AtomType(atom_name, Mass=1., Charge = charge)
        sim_atom_types[atom_name] = atom_type
    print('atom types: {}'.format(sim_atom_types))

    # 2b) define molecules
    for mol_name in topology.molecule_names:
        if sys_params.num_mols[mol_name]>0:
            atoms_in_mol = []
            #get list of atom types
            for atom_name in topology.molecule_defs[mol_name]:
                atoms_in_mol.append( sim_atom_types[atom_name] )
            #initialize/make the moltype
            mol_type = sim.chem.MolType(mol_name, atoms_in_mol)
            #add bonds
            if topology.bonding is not None:
                if mol_name in topology.bonding:
                    for bond in topology.bonding[mol_name]:
                        mol_type.Bond(bond[0],bond[1])
            #append
            sim_mol_types.append(mol_type)
    print('molecule types: {}'.format(sim_mol_types))

    # 2c) create 'World' and system
    World = sim.chem.World(sim_mol_types, Dim=3, Units=sys_params.units)
    sys_name = special_options.sys_name if (special_options and special_options.sys_name) else sys_params.sys_name
    Sys = sim.system.System(World, Name = sys_name)
    Sys.BoxL = sys_params.box
    print('Setting box: {}'.format(Sys.BoxL))

    for i, mol_type in enumerate(sim_mol_types):
        num_mol = sys_params.num_mols[mol_type.Name]

        print("Adding {} {} molecules to system".format(num_mol, mol_type.Name))
        for j in range(num_mol):
            Sys += mol_type.New()

    # 3) create forcefield
    print('... Initializing Force Field ...')
    force_field = ff.create_forcefield(Sys,topology,interactions, special_options)
    Sys.ForceField.extend(force_field)
    print('Force Field Summary: {}'.format(Sys.ForceField))

    if load:
        load_system(Sys,ff_file=sys_params['forcefield_file'])
    """
        # set up initial parameters
        if sys_params['forcefield_file']: 
            #with open(sys_params['forcefield_file'], 'r') as of: s = of.read()
            #Sys.ForceField.SetParamString(s)      
            ff.set_param_from_file(Sys.ForceField, sys_params['forcefield_file'])

        # set up histograms
        for P in Sys.ForceField:
            P.Arg.SetupHist(NBin = 10000, ReportNBin = 100)
        # lock and load
        Sys.Load()
    """
   
    # 4) Other options
    print('... Setting other options, e.g. temperature and integration ...')
    Sys.TempSet = sys_params['temp']
    Sys.PresSet = sys_params['pressure']
    Sys.PresAx = sys_params['barostat_options']['pressure_axis']
    Sys.BarostatOptions['tension'] = sys_params['barostat_options']['tension']
    
    Int = Sys.Int
    Int.Method = Int.Methods.VVIntegrate        
    Int.Method.TimeStep = sys_params['dt']
    if Sys.PresSet <= 0.: #NVT
        Int.Method.Thermostat = Int.Method.ThermostatLangevin
        Int.Method.LangevinGamma = sys_params['langevin_gamma']
    else: #NPT
        Int.Method.Thermostat = Int.Method.ThermostatNoseHoover
        Int.Method.Barostat = Int.Method.BarostatMonteCarlo  

    return Sys

def load_system(_Sys, ff_file=None):
    # set up initial parameters
    if ff_file is not None: 
        #with open(sys_params['forcefield_file'], 'r') as of: s = of.read()
        #Sys.ForceField.SetParamString(s)      
        ff.set_param_from_file(_Sys.ForceField, ff_file)
   
    print('--- Force Field Summary: {} ---'.format(_Sys.ForceField))
    print('{}'.format(_Sys.ForceField.ParamString()))

    # set up histograms
    for P in _Sys.ForceField:
        P.Arg.SetupHist(NBin = 10000, ReportNBin = 100)
    # lock and load
    _Sys.Load()


