# Helper functions for setting up a sim object from my custom "universal" data structures
# (C) Kevin Shen, 2021
#
# TODO:
# 1) in future, along with forcefield.py, consider breaking out each potential into its own module to gather ff parsing, definition, and exports all in one place
# 2) writing from Sys.ForceField back out to my ff_dict format
#     only have simple implementation right now. need more work to get it to write a new forcefield object, de novo! tackle later once I have force field shorthand export as well.
#   to export force field fully, would need:
#   1) function to determine species AtomTypes
#   2) know what parameters to extract from the potential. May end up needing to define function_specific functions, like we did for setup()
#   3) if section doesn't exist, make it. fftype:params_sim:listofpotentials
#      can make it export in shorthand, just need list of strings, and then two join commands!
#
#
# Important Usage Notes:
# 1) requires dictionary file format to be in my specified format, with naming conventions for particular fields.
#    i.e. relies on forcefield.ForceField() to do the parsing and making sure the important fields are in there
#    i.e. having a params_sim section, naming convention for parameters and potential types, etc.
#    in future, can think about how to better modularize all of these naming conventions
#
#
import sys,os
#sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from context import sim
import sim
import re
import copy
import topologify,forcefield
from collections import OrderedDict

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
  sys_name = options['sys_name'] if (options and 'sys_name' in options and options['sys_name']) else 'Sys'
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
  options used:
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
    - neutralize
  '''
  print('===== Creating ForceField =====')
  if isinstance(ff_dict, forcefield.ForceField):
    ff_dict = ff_dict.processed_file

  sim_atom_types = { atom.Name:atom for atom in Sys.World.AtomTypes }
  forcefields = []

  for ff_type,setup_func in setup_dict.items():
    print('---> Setting up ff type {}'.format(ff_type))
    if ff_type == 'coulomb_smeared' and 'neutralize' in options and options['neutralize'] == True:
      pass
    elif ff_type in ff_dict:
      forcefields.extend( setup_func( Sys, sim_atom_types, ff_dict[ff_type] ) )
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
  param_names = ['Dist0','FConst']
  addtl_dict = {'Fixed':True} #just default. fixables will be updated by ff_dict.
  ffs = setup_ff_base( Sys, sim_atom_types, ff_dict, ff_type, param_names, addtl_dict )
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
  param_names = ['Cut','B','Kappa','Dist0','Sigma','Epsilon']
  addtl_dict = {'Fixed':True,'Shift':False} #just default. fixables will be updated by ff_dict.
  ffs = setup_ff_base( Sys, sim_atom_types, ff_dict, ff_type, param_names, addtl_dict )
  #my default is to set B.Min to negative 100
  for ff in ffs:
    print('setting B.Min to -100.0 for {}'.format(ff.Label))
    ff.B.Min = -100.0
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
  param_names = ['Cut','Coef','BornA']
  addtl_dict = {'FixedCoef':True, 'FixedBornA':True} #just default. fixables will be updated by ff_dict.
  ffs = setup_ff_base( Sys, sim_atom_types, ff_dict, ff_type, param_names, addtl_dict )

  # No matter what, need to add Ewald
  ff_defaults = ff_dict['defaults_sim']
  p = sim.potential.Ewald( Sys, ExcludeBondOrd=int(ff_defaults['ExcludeBondOrd']), 
      Cut=ff_defaults['Cut']['val'], Shift=ff_defaults['Shift'], Coef=ff_defaults['Coef']['val'], FixedCoef=ff_defaults['Coef']['fixed'], Label='ewald' )
  ffs.append(p)

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
  param_names = ['UConst', 'NPeriods', 'PlaneAxis', 'PlaneLoc']
  addtl_dict = {'Fixed':True} #just default. fixables will be updated by ff_dict.
  ffs = setup_ff_base( Sys, sim_atom_types, ff_dict, ff_type, param_names, addtl_dict )
  return ffs   

def setup_ff_base( Sys, sim_atom_types, ff_dict, ff_type, param_names, addtl_dict ):
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
  param_names : list
    list of argument names that can be retrieved in ff_dict['param_sim']
  addtl_dict : dict
    additional arguments special to the potential type that need to be set, but are not in ff_dict
    

  Notes
  -----
  relies on the dictionaries being fed in having all the fields that are required.
  IMPORTANT naming conventions/requirements: 
    requires is_bonded_dict to know whether ff_type is bonded or not
    requires sim_potential_dict to know what sim.potential to call with given ff_type
    requires param_names in ff_dict to be same as in sim definition! otherwise, will need a naming mapping dict to know how to parse ff_dict. can maybe implement later...
  '''

  ffs = []
  # Begin processing
  for ff_entry in ff_dict['params_sim']:
    # first check atom types are valid
    if atomtypes_in_system( ff_entry['species'], sim_atom_types ):
      print('adding {} for species {}'.format(ff_type, ff_entry['species']))
      atypes = [ [sim_atom_types[aname] for aname in species_set] for species_set in ff_entry['species'] ]
      if len(atypes) == 1:
        f = sim.atomselect.Filter( atypes )
      elif len(atypes) == 2:
        bonded = is_bonded_dict[ff_type]
        if bonded: f = sim.atomselect.PolyFilter( atypes, Bonded=bonded )
        else: f = sim.atomselect.PolyFilter( atypes ) #Bonded=None

      #define
      param_dict = { param_name:ff_entry[param_name]['val'] for param_name in param_names }
      for k,v in addtl_dict.items(): #add the additional default arguments to param_dict
        if k in ff_entry:
          param_dict[k] = ff_entry[k]
        else:
          param_dict[k] = v
      p = sim_potential_dict[ff_type]( Sys, Filter=f, Label = ff_entry['name'], **param_dict)
      print(param_dict)
      print(ff_entry['name'])
      #fixable_dict= { name:ii for ii,name in enumerate(p.Param.Names) }
      #for param_name in fixable_dict:
      #  p.Param.Fixed[fixable_dict[param_name]] = ff_entry[param_name]['fixed']
      for param_name in p.Param.Names:
        p.__getattr__(param_name).Fixed = ff_entry[param_name]['fixed']

      #collect
      ffs.append(p)
    else:
      raise ValueError('Some atom types required by {} potential {} not present, skipping'.format( ff_type, ff_entry ) )
  return ffs   


setup_dict = OrderedDict( [('bond_harmonic',setup_bond_harmonic), ('pair_ljg',setup_pair_ljg),
              ('coulomb_smeared',setup_coulomb_smeared), ('external_sin',setup_external_sin)] )
sim_potential_dict = { 'bond_harmonic':sim.potential.Bond, 
    'pair_ljg':sim.potential.LJGaussian,
    'external_sin':sim.potential.ExternalSinusoid, 
    'coulomb_smeared':sim.potential.SmearedCoulombEwCorr }
is_bonded_dict = { 'bond_harmonic': True, 'external_sin':False, 'coulomb_smeared':False, 'pair_ljg':False }


def save_params_to_dict( force_field, out_dict ):
  ''' save parameters to ff_dict
  Parameters
  ----------
  force_field: sim.potential.ForceField list
  out_dict : dict-like
    Should be the root-level, i.e. full system force field definition. Will save parameters to this dict.

  Notes
  -----
  Will save using sim's variable naming convention. I.e. the output dictionary should have parameters in that format as well.
  Version 0: Assume that the out_dict already has the entire structure set up already, i.e. all fields already exist
  Version 1: Create and format the parameters appropriately in out_dict if they don't exist. Note that some of the common parameters, like cutoff, etc. won't be included in the first pass.
  '''
  ffdict = {p.Name:p for p in force_field}
  fftype_map = {'bond':'bond_harmonic', 'ljg':'pair_ljg', 
        'smearedcoulombEwCorr':'coulomb_smeared', 
        'external_sinusoid':'external_sin', 'ewald':'skip'}

  # For efficiency, collect potentials by type
  ff_by_type = {}
  for ff_sim in force_field:
    ff_type_sim = ff_sim.Names[0]
    if ff_type_sim not in ff_by_type:
      ff_by_type[ ff_type_sim ] = []
    ff_by_type[ ff_type_sim ].append(ff_sim)

  # Now go through and set the potentials
  for ff_type_sim,ffs in ff_by_type.items():
    ff_type_out = fftype_map[ ff_type_sim ]
    if ff_type_out in ['skip']:
      continue
    output_fflist = out_dict[ff_type_out]['params_sim'] 
    ff_by_index_in_outdict = { ff['name']:ii for ii,ff in enumerate(output_fflist)}
    print('---')
    print(output_fflist)
    print(ff_by_index_in_outdict)
    for ff_sim in ffs:
      ff_out = output_fflist[ ff_by_index_in_outdict[ff_sim.Label] ]
      for param_name in ff_sim.Param.Names:
        ff_out[param_name]['val'] = ff_sim.__getattr__(param_name)[0]


# ===== misc. helpers for working with sim =====
def set_params_from_file(force_field, ff_file):
  ''' Set sim ForceField list from a file
  Parameters
  ----------
  force_field: sim.potential.ForceField list
  ff_file : str
    ff_file to load
  '''
  with open(ff_file, 'r') as of: s = of.read()
  print('... Setting FF with file {} ...\n{}'.format(ff_file,s))
  force_field.SetParamString(s)      

def update_specific_params(force_field, params):
  """Update specific terms of an existing Sys.ForceField
  Parameters
  ----------
  force_field: sim.potential.ForceField list
  params: list or tuple

  Notes
  -----
  Expecting format is ('Name',bool) to set default for an entire potential
  ('Name',(parameter,value,bool)) to set a specific parameter in a specific forcefield
  
  Also, this doesn't check for conflicts in the input params! I.e. this will be first in first out, later settings will override earlier settings
  """
  # first create look-up dict of forcefield
  ffdict = {p.Name:p for p in force_field}

  def set_parameter(param):
    potential_name = param[0]
    if isinstance(param[1],bool):
      print('Setting {} fixed status to all {}'.format(potential_name,param[1]))
      ffdict[potential_name].Fixed[:] = param[1]
    elif isinstance(param[1],(list,tuple)):
      print('Setting {}, {}'.format(potential_name,param[1]))
      param_name = param[1][0]
      for entry in param[1][1:]:
        if isinstance(entry,bool):
          ffdict[potential_name].__getattr__(param_name).Fixed = entry
        elif isinstance(entry,float):
          ffdict[potential_name].__getattr__(param_name)[0] = entry 
    else:
      raise ValueError("given param {} option not understood".format(param))

  if isinstance(params,(list,tuple)):
    if not isinstance(params[0],(list,tuple)):
      param = params
      set_parameter(param)
    else:
      print('detected list of tuples')
      for param in params: 
        set_parameter(param)       
  else: 
    raise ValueError('params should be entered as tuple/list of tuples and lists')

def get_free_params(force_field):
  ''' Get list of currently free parameters 
  Parameters
  ----------
  force_field: sim.potential.ForceField list

  Returns
  -------
  free_params: list
    list of (potential.Name, param_name) tuples for easy access later
  '''
  free_params = []
  for potential in force_field:
    for ip,param_name in enumerate(potential.Param.Names):
      fixed = potential.Fixed[ip]
      if not fixed:
        free_params.append( (potential.Name,param_name) )
  return free_params

def toggle_params_fixed(force_field,params,fixed):
  ''' toggle the fixed state of specified parameters to a specified boolean value
  Parameters
  ----------
  force_field : sim.potential.ForceField list
  params : list
    list of (potential.Name, param_name) tuples for easy access later
  fixed: boolean
    the boolean value to toggle to

  Notes
  -----
  In the future can consider the default toggle state to switch to opposite value, instead of blanket setting the value of everything in `params` to `fixed`
  '''
  ffdict = {p.Name:p for p in force_field}
  for param in params:
    potential_name = param[0]
    parameter_name = param[1]
    ffdict[potential_name].__getattr__(parameter_name).Fixed = fixed

def load_system( _Sys, ff_file = None ):
  ''' loads Sys object, i.e. compiles Fortran and locks it in. 
  Parameters
  ----------
  _Sys : sim.system
  ff_file : str
    sim-format force field file to load

  Notes
  -----
  Also sets up histogram. Can consider exposing those default settings somewhere else.
  As an idiom/shorthand, allowing specifying ff_file with this function, s.t. can load and get a system set up for optimizing in one go.
  '''
  # set up initial parameters
  if ff_file is not None: 
      set_params_from_file(_Sys.ForceField, ff_file)
 
  print('--- Force Field Summary: {} ---'.format(_Sys.ForceField))
  print('{}'.format(_Sys.ForceField.ParamString()))

  # set up histograms
  for P in _Sys.ForceField:
      P.Arg.SetupHist(NBin = 10000, ReportNBin = 100)

  # lock and load
  _Sys.Load()

#usually, for saving, just use Opt's wrapper for Sys.ForceField.ParamString()
#def save_system(_Sys):
#  ''' save force field for system '''



# ===== Notes =====
# Set up optimizer


