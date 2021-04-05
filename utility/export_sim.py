# Helper functions for setting up a sim object from my custom "universal" data structures
# (C) Kevin Shen, 2021
#
# TODO:
# 1) in future, along with forcefield.py, consider breaking otu each potential into its own module to gather ff parsing, definition, and exports all in one place
#
import sys
sys.path.insert(1,'../')
import sim
import re
import copy
import topologify,forcefield

## ===== parsers and utilities =====
def isiterable( maybe_iterable ):
  '''https://stackoverflow.com/questions/1952464/in-python-how-do-i-determine-if-an-object-is-iterable'''
  try: 
    iter( maybe_iterable )
    return True
  except:
    return False

def atomtypes_in_system(atype_names, sim_atom_types):
  '''
  atype_names: list of lists
  sim_atom_types: dict,list,tuple
    {atomtypename: simatomtype object} or simply list/tuple of atom type names

  Notes
  -----
  atype_names is expected to be of the format:
    [ [typename1, typename2,...] , [ typename1,...] ] for 2-body
    [ [typename1, typename2,...] ] for 1-body, note is still a list of list.
  '''
  processed = []
  for element in atype_names:
    if isinstance(element,list) or isinstance(element,tuple):
      processed.append( atomtypes_in_system(element,sim_atom_types) )
    else:
      processed.append(element in sim_atom_types)
  return all(processed)




## ===== system =====
# Set up topology (or chemistry/world/system in sim's language)
def create_system( mytop, options=None ):
  '''
  Parameters
  ----------
  mytop : topologify.Topology
    or dict of the form in topologify.Topology.processed_file (i.e. no shorthand)
  options : dict
    fields used: neutralize, sys_name,

  Returns
  -------
  sys : sim Sys
  '''
  print('===== Creating System (Topology) =====')
  units = sim.units.DimensionlessUnits
  if isinstance(mytop,topologify.Topology):
    mytop = mytop.processed_file

  # 1) Define atom and mol types
  sim_atom_types = {}
  sim_mol_types_list = []
  sim_mol_types_dict = {}

  print('---> Defining atom types')
  for atom in mytop['bead_types']:
    if 'name' in atom: #otherwise is a default field specification
      atom_name = atom['name']
      mass = atom['mass']
      charge = 0.0 if (options and 'neutralize' in options and options['neutralize']) else atom['charge']

      atom_type = sim.chem.AtomType(atom_name, Mass=mass, Charge = charge)
      sim_atom_types[atom_name] = atom_type 
  print('atom types: {}\n'.format(sim_atom_types))

  print('---> Defining molecule types')
  for mol in mytop['mol_types']:
    atoms_in_mol = [ sim_atom_types[atom_name] for atom_name in mol['beads'] ]

    mol_type = sim.chem.MolType( mol['name'], atoms_in_mol )
    for bond in mol['bonds']: #stored as intra-chain site index
      mol_type.Bond( bond[0],bond[1] )

    sim_mol_types_list.append(mol_type)
    sim_mol_types_dict[ mol['name'] ] = mol_type
  print('molecule types: {}\n'.format(sim_mol_types_list))
    
  # 2) Create 'World' and system 
  print('---> Creating World and Sys')
  World = sim.chem.World(sim_mol_types_list, Dim=3, Units=units)
  sys_name = options['sys_name'] if (options and 'sys_name' in otpions and options['sys_name']) else 'Sys'
  Sys = sim.system.System(World, Name = sys_name)

  for im, entry in enumerate( mytop['system'] ):
    mol_name = entry[0]
    num_mol = entry[1]
    mol_type = sim_mol_types_dict[ mol_name ]

    print("Adding {} {} molecules to system".format(num_mol, mol_name))
    for j in range(num_mol):
      Sys += mol_type.New()
  print('')

  # Final stuff
  return Sys

def set_system_settings( Sys, options ):
  ''' set up settings (ensemble, integrator, etc.)
  Notes
  -----
  Preliminary implementation: require all the fields to be present to have well-defined system
  '''
  #box:
  Sys.BoxL = options['box'] #should be 3-element listable object
  print('Setting box: {}'.format(Sys.BoxL))

  #ensemble and integrator:
  Sys.TempSet = options['temp']
  if 'barostat' in options:
    if 'pressure' in options['barostat']:
      Sys.PresSet = options['barostat']['pressure']
    if 'pressure_axis' in options['barostat']:
      Sys.PresAx = options['barostat']['pressure_axis']
    if 'tension' in options['barostat']:
      Sys.BarostatOptions['tension'] = options['barostat']['tension']
  
  Int = Sys.Int
  Int.Method = Int.Methods.VVIntegrate
  Int.Method.TimeStep = options['integrator']['dt']
  if 'langevin_gamma' in options['integrator']:
    Int.Method.LangevinGamma = options['integrator']['langevin_gamma']
  if Sys.PresSet <= 0. :#NVT
    Int.Method.Thermostat = Int.Method.ThermostatLangevin
    Int.Method.Barostat = Int.Method.BarostatNone
  else:
    Int.Method.Thermostat = Int.Method.ThermostatNoseHoover
    Int.Method.Barostat = Int.Method.BarostatMonteCarlo  

  return Sys

# ===== Set up force field =====
def create_forcefield( Sys, mytop, ff_dict, options={} ):
  '''
  Parameters
  ----------
  Sys : sim System object
  mytop : topologify.Topology()
  ffdict : dictionary-like
    should be processed, no shorthand

  Notes
  -----
  options: special options flags

  '''
  print('===== Creating ForceField =====')
  if isinstance(ff_dict, forcefield.ForceField):
    ff_dict = ff_dict.processed_file
  setup_dict = {'bond_harmonic':setup_bond_harmonic, 'pair_ljg':setup_pair_ljg,
                'coulomb_smeared':setup_coulomb_smeared, 'external_sin':setup_external_sin}
  sim_atom_types = { atom.Name:atom for atom in Sys.World.AtomTypes }
  forcefields = []

  for ff_type,setup in setup_dict.items():
    print('---> Setting up ff type {}'.format(ff_type))
    forcefields.append( setup( Sys, sim_atom_types, ff_dict[ff_type] ) )
    print('')

  return forcefields


def setup_bond_harmonic( Sys, sim_atom_types, ff_dict ):
  ''' set up harmonic bonds
  Parameters
  ----------
  Sys: sim System object
  sim_atom_types : dict
    {atom_type_name: sim_atom_type_object, ...}
  ff_dict: dictionary
    the ff.processed_file['bond_harmonic'] section/format
  '''
  ff_type = 'bond_harmonic'
  ffs = []
  # Begin processing
  for ff_entry in ff_dict['params_sim']:
    # first check atom types are valid
    if atomtypes_in_system( ff_entry['species'], sim_atom_types ):
      print('adding {} for species {}'.format(ff_type, ff_entry['species']))
      atypes0 = [ sim_atom_types[aname] for aname in ff_entry['species'][0] ]
      atypes1 = [ sim_atom_types[aname] for aname in ff_entry['species'][1] ]
      f = sim.atomselect.PolyFilter( [atypes0,atypes1], Bonded=True )

      #define
      p = sim.potential.Bond(Sys, Filter=f, Fixed=True, Dist0=ff_entry['Dist0']['val'], FConst=ff_entry['FConst']['val'], Label=ff_entry['name'])
      p.Param.Dist0.Fixed = ff_entry['Dist0']['fixed']
      p.Param.FConst.Fixed = ff_entry['FConst']['fixed']

      #collect
      ffs.append(p)
    else:
      raise ValueError('Some atom types required by {} potential {} not present, skipping'.format( ff_type, ff_entry ) )
  return ffs 

 
def setup_pair_ljg( Sys, sim_atom_types, ff_dict ):
  ''' set up LJG interaction
  Parameters
  ----------
  Sys: sim System object
  sim_atom_types : dict
    {atom_type_name: sim_atom_type_object, ...}
  ff_dict: dictionary
    the ff.processed_file['pair_ljg'] section/format
  '''
  ff_type = 'pair_ljg'
  ffs = []
  # Begin processing
  for ff_entry in ff_dict['params_sim']:
    # first check atom types are valid
    if atomtypes_in_system( ff_entry['species'], sim_atom_types ):
      print('adding {} for species {}'.format(ff_type, ff_entry['species']))
      atypes0 = [ sim_atom_types[aname] for aname in ff_entry['species'][0] ]
      atypes1 = [ sim_atom_types[aname] for aname in ff_entry['species'][1] ]
      f = sim.atomselect.PolyFilter( [atypes0,atypes1], Bonded=False )

      #define
      p = sim.potential.LJGaussian( Sys, Filter=f, Cut=ff_entry['Cut']['val'], Fixed=True,
              B=ff_entry['B']['val'], Kappa=ff_entry['Kappa']['val'], Dist0=ff_entry['Dist0']['val'], 
              Sigma=ff_entry['Sigma']['val'], Epsilon=ff_entry['Epsilon']['val'], Label=ff_entry['name'] )

      p.Param.B.Fixed = ff_entry['B']['fixed']
      p.Param.B.Min = -100. #default, don't restrict positive
      p.Param.Kappa.Fixed = ff_entry['Kappa']['fixed']
      p.Param.Dist0.Fixed = ff_entry['Dist0']['fixed']
      p.Param.Sigma.Fixed = ff_entry['Sigma']['fixed']
      p.Param.Epsilon.Fixed = ff_entry['Epsilon']['fixed']

      #collect
      ffs.append(p)
    else:
      raise ValueError('Some atom types required by {} potential {} not present, skipping'.format( ff_type, ff_entry ) )
  return ffs 

def setup_coulomb_smeared( Sys, sim_atom_types, ff_dict ):
  ''' set up smeared coulomb potential
  Parameters
  ----------
  Sys: sim System object
  sim_atom_types : dict
    {atom_type_name: sim_atom_type_object, ...}
  ff_dict: dictionary
    the ff.processed_file['coulomb_smeared'] section/format

  Notes
  -----
  Note Ewald needs periodic boundary condition, *requires* box size to be set already!
  '''
  ff_type = 'coulomb_smeared'
  ffs = []

  # No matter what, need to add Ewald
  ff_defaults = ff_dict['defaults_sim']
  p = sim.potential.Ewald( Sys, ExcludeBondOrd=int(ff_defaults['ExcludeBondOrd']), 
      Cut=ff_defaults['Cut']['val'], Shift=ff_defaults['Shift'], Coef=ff_defaults['Coef']['val'], FixedCoef=ff_defaults['Coef']['fixed'], Label='ewald' )
  ffs.append(p)

  # Begin processing
  for ff_entry in ff_dict['params_sim']:
    # first check atom types are valid
    if atomtypes_in_system( ff_entry['species'], sim_atom_types ):
      print('adding {} for species {}'.format(ff_type, ff_entry['species']))
      atypes0 = [ sim_atom_types[aname] for aname in ff_entry['species'][0] ]
      atypes1 = [ sim_atom_types[aname] for aname in ff_entry['species'][1] ]
      f = sim.atomselect.PolyFilter( [atypes0,atypes1], Bonded=False )

      #define
      p = sim.potential.SmearedCoulombEwCorr( Sys, Filter=f, 
          Cut=ff_entry['Cut']['val'], Shift=ff_entry['Cut']['val'], 
          Coef=ff_entry['Coef']['val'], FixedCoef=ff_entry['Coef']['fixed'], 
          BornA=ff_entry['BornA']['val'], FixedBornA=ff_entry['BornA']['fixed'], 
          Label=ff_entry['name'])

      #collect
      ffs.append(p)
    else:
      raise ValueError('Some atom types required by {} potential {} not present, skipping'.format( ff_type, ff_entry ) )
  return ffs   


def setup_external_sin( Sys, sim_atom_types, ff_dict ):
  ''' set up external sin potential
  Parameters
  ----------
  Sys: sim System object
  sim_atom_types : dict
    {atom_type_name: sim_atom_type_object, ...}
  ff_dict: dictionary
    the ff.processed_file['external_sin'] section/format
  '''
  ff_type = 'external_sin'
  ffs = []
  # Begin processing
  for ff_entry in ff_dict['params_sim']:
    # first check atom types are valid
    if atomtypes_in_system( ff_entry['species'], sim_atom_types ):
      print('adding {} for species {}'.format(ff_type, ff_entry['species']))
      atypes0 = [ sim_atom_types[aname] for aname in ff_entry['species'][0] ]
      f = sim.atomselect.Filter( atypes0 )

      #define
      paramdict = { 'UConst':ff_entry['UConst']['val'], 'NPeriods':ff_entry['NPeriods']['val'],
          'PlaneAxis':ff_entry['PlaneAxis']['val'], 'PlaneLoc':ff_entry['PlaneLoc']['val'], 
          'Fixed': True
          }
      p = sim.potential.ExternalSinusoid( Sys, Filter=f, Label = ff_entry['name'], **paramdict)
      for pname in paramindex:
        p.Param.Fixed[paramindex[pname]] = ff_entry[pname]['fixed']

      '''
      p = sim.potential.ExternalSinusoid(Sys, Filter=f, Fixed=True, 
           UConst=ff_entry['UConst']['val'], NPeriods=ff_entry['NPeriods']['val'],
           PlaneAxis=ff_entry['PlaneAxis']['val'], PlaneLoc=ff_entry['PlaneLoc']['val'], 
           Label=ff_entry['name'])
      p.Param.UConst.Fixed = ff_entry['UConst']['fixed']
      p.Param.NPeriods.Fixed = ff_entry['NPeriods']['fixed']
      '''

      #collect
      ffs.append(p)
    else:
      raise ValueError('Some atom types required by {} potential {} not present, skipping'.format( ff_type, ff_entry ) )
  return ffs   


def setup_generic( Sys, sim_atom_types, ff_dict, ff_type, param_dict ):
  ''' set up generic potetnails, currently for 1 or 2-species defined interactions
  Parameters
  ----------
  Sys: sim System object
  sim_atom_types : dict
    {atom_type_name: sim_atom_type_object, ...}
  ff_dict: dictionary
    the ff.processed_file[ ff_type_name ] section/format
  ff_type : str
    name of the type of potential, for documentation purposes only
  param_dict : dict
    dict of arguments that are to be fed into the potential constructor, { param_name:param_value }

  Notes
  -----
  relies on the dictionaries being fed in having all the fields that are required.
  '''
  ffs = []
  # Begin processing
  for ff_entry in ff_dict['params_sim']:
    # first check atom types are valid
    if atomtypes_in_system( ff_entry['species'], sim_atom_types ):
      print('adding {} for species {}'.format(ff_type, ff_entry['species']))
      atypes0 = [ sim_atom_types[aname] for aname in ff_entry['species'][0] ]
      f = sim.atomselect.Filter( atypes0 )

      #define
      p = sim.potential.ExternalSinusoid( Sys, Filter=f, Label = ff_entry['name'], **paramdict)
      fixable_dict= { name:ii for ii,name in enumerate(p.Param.Names) }
      for param_name in fixable_dict:
        p.Param.Fixed[paramindex[param_name]] = ff_entry[param_name]['fixed']

      #collect
      ffs.append(p)
    else:
      raise ValueError('Some atom types required by {} potential {} not present, skipping'.format( ff_type, ff_entry ) )
  return ffs   

# ===== misc. helpers =====
def load_system( Sys, ff_file ):
  ''' loads Sys object, i.e. compiles Fortran and locks it in. '''
  pass



# ===== Notes =====
# Sys.ForceField.extend(force_field)
# load_system

# Other updates and settings
#Still required:
'''
Sys.BoxL = sys_params.box
print('Setting box: {}'.format(Sys.BoxL))
'''

# Set up optimizer


# Set up other stuff
