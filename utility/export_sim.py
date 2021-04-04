# Helper functions for setting up a sim object from my custom "universal" data structures
# (C) Kevin Shen, 2021
#
import sys
sys.path.insert(1,'../')
import sim
import re
import copy

## parsers

## system
# Set up topology (or chemistry/world/system in sim's language)
def create_system( mytop, options=None ):
  '''
  Parameters
  ----------
  mytop : topologify Topology
  options : dict
    fields used: neutralize, sys_name,

  Returns
  -------
  sys : sim Sys
  '''
  units = sim.units.DimensionlessUnits

  # 1) Define atom and mol types
  sim_atom_types = {}
  sim_mol_types_list = []
  sim_mol_types_dict

  for atom in mytop.processed_file['beads']:
    atom_name = atom['name']
    mass = atom['mass']
    charge = 0.0 if (options and options['neutralize']) else atom['charge']

    atom_type = sim.chem.AtomType(atom_name, Mass=mass, Charge = charge)
    sim_atom_types[atom_name] = atom_type 
  print('atom types: {}'.format(sim_atom_types))

  for mol in mytop.processed_file['mol_type']:
    atoms_in_mol = [ sim_atom_types[atom_name] for atom_name in mol['beads'] ]

    mol_type = sim.chem.MolType( mol['name'], atoms_in_mol )
    for bond in mol['bonds']: #stored as intra-chain site index
      mol_type.Bond( bond[0],bond[1] )

    sim_mol_types_list.append(mol_type)
    sim_mol_types_dict[ mol['name'] ] = mol_type
  print('molecule types: {}'.format(sim_mol_types_list))
    
  # 2) Create 'World' and system 
  World = sim.chem.World(sim_mol_types_list, Dim=3, Units=units)
  sys_name = options['sys_name'] if (options and options['sys_name']) else 'Sys'
  Sys = sim.system.System(World, Name = sys_name)

  for im, entry in enumerate( mytop.processed_file['system'] ):
    mol_name = entry[0]
    num_mol = entry[1]
    mol_type = sim_mol_types_dict[ mol_name ]

    print("Adding {} {} molecules to system".format(num_mol, mol_name))
    for j in range(num_mol):
      Sys += mol_type.New()
 
  # Final stuff
  return Sys


# ===== Set up force field =====
def create_forcefield( Sys, mytop, interactions, options=None ):
  return 0 

def parse_ljg( ljg_dict ):
  '''
  Take in shorthand, and give processed, fully verbose form.
  '''
  for line_entry in ljg_dict['params_sim']:
    return 0
   


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
