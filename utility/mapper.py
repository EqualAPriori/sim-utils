# (C) Kevin Shen, 2021, kevinshen@ucsb.edu
# Utility functions for doing mappings
# 
# Main function types:
# 1) I/O
# 2) utilities for generating the index mappings
# 3) utililities for actually doing the mapping of a configuration/trajectory
#
import numpy as np
import mdtraj
import ruamel.yaml as YAML
import os

# ===== File I/O =====
def create_yaml():
    '''
    Notes
    -----
    https://stackoverflow.com/questions/49669236/ruamel-yaml-bad-dump
    '''
    yaml = YAML.YAML()
    yaml.explicit_start = True
    yaml.default_flow_style = None 
    yaml.encoding = "utf-8"     # default when using YAML() or YAML(typ="rt")
    yaml.allow_unicode = True   # always default in the new API
    yaml.errors = "strict"
    return yaml
yaml = create_yaml()

def save(prefix,aa_indices_in_cg,cg_site_of_aa,shorthand=None,traj_mapped=None, traj_mapping=None):
    mappingfile = prefix + '_mapping.yaml'
    if traj_mapped is not None:
        pdbfile = prefix + '_mapped.pdb'
        traj_mapped.save(pdbfile)
    if traj_mapping is not None:
        pdbfile_mapping = prefix + '_mapping.pdb'
        traj_mapping.save(pdbfile_mapping)
    else: pdbfile_mapping = None
    d = {
            'pdbfile_mapped': pdbfile,
            'pdbfile_mapping': pdbfile_mapping,
            'aa_indices_in_cg': aa_indices_in_cg,
            'cg_site_of_aa': cg_site_of_aa,
            'shorthand': shorthand
        }

    with open(mappingfile,'w') as f:
        f.write('# mapping specifications for molecule {}\n'.format(prefix))
        yaml.dump(d,f)


# ===== Mapping Utilities =====
def generate_pdb_mapping(traj,mapping):
    '''
    take in a cg_site_of_aa mapping array, and generate a new mapping topology for the chain, annotated with cgbead site
    Parameters
    ----------
    traj : mdtraj trajectory
    mapping : array
        should be cg_site_of_aa array format, of a single chain
    '''
    mapping_top = mdtraj.Topology()
    ch_previous = None
    sites = []
    site_to_res_map = {}
    for ia,atom in enumerate(traj.top.atoms):
        cg_beadname, cg_identifier = mapping[ia]
        if ia == 0:
            #if atom.chain != ch_previous: #techincally I envision this only working for single chain first
            ch_previous = mapping_top.add_chain()
        if atom.residue.chain != traj.top.atom(0).residue.chain:
            raise ValueError('Multiple chains detected, not implemented yet')
        if cg_identifier not in sites: #need to add a new residue
            sites.append(cg_identifier)
            res = mapping_top.add_residue(cg_beadname,ch_previous)
            site_to_res_map[cg_identifier] = res
        
        res = site_to_res_map[cg_identifier]
        new_atom = mapping_top.add_atom(atom.name, atom.element, res)

    for bond in traj.top.bonds:
        a1,a2 = mapping_top.atom(bond[0].index), mapping_top.atom(bond[1].index)
        mapping_top.add_bond( a1,a2 )

    mapping_traj = mdtraj.Trajectory( traj.xyz, mapping_top )
    return mapping_traj
 
def generate_single_mapping_shorthand(mode,shorthand):
    '''
    A helper function to generate requisite data structure for generate_single_mapping() from a more compact notation.

    Parameters
    ----------
    shorthand : str or list
    mode = option

    Returns
    -------
    aa_indices_in_cg : list of lists
        aa_indices_in_cg[i] returns list of aa indices in cg bead i
    cg_site_of_aa : list
        cg_site_of_aa[i] returns the cg site that atom i maps to, for us in generate_single_mapping

    Notes
    -----
    Shorthand versions:
    1) assuming AA atoms are in contiguous blocks for mapping
       [ naa1, naa2, ... ] each entry is # AA atoms in the next CG bead
    2) ... smiles? some other graph format? 
    3) I guess can technically wrap this up into generate_single_mapping
    '''
    if mode.lower() in ['simple','contiguous','short']:
        # expect data to be of the form [(cg_bead_name, n_aa_in_cg)]
        #cg_site_of_aa = [ ii for ii,num_in_cgbead in enumerate(shorthand) for jj in range(num_in_cgbead) ]
        nbead_so_far = 0
        aa_indices_in_cg = []
        cg_site_of_aa = []
        for cgindex,entry in enumerate(shorthand):
            if (isinstance(entry,list) or isinstance(entry,tuple)):
                assert len(entry) == 2, "mapping entry for simple shorthand should be an integer or a 2-element tuple"
                cg_bead_name,num_in_cgbead = entry
            elif isinstance(entry,int):
                cg_bead_name,num_in_cgbead = 'DMY',entry 
            else:
                raise ValueError('Unrecognized mapping entry for simple shorthand: {}'.format(entry))
            indices = np.arange(nbead_so_far,nbead_so_far + num_in_cgbead).tolist()
            aa_indices_in_cg.append( (cg_bead_name,indices) ) 
            nbead_so_far += num_in_cgbead
            for ii in range(num_in_cgbead):
                cg_site_of_aa.append( (cg_bead_name,cgindex) )
        #cg_site_of_aa = []
        #aa_indices_in_cg = []
        #for ii,cgbead in enumerate(shorthand):
        #    aa_indices_in_cg.append( [ii]*cgbead )
        #    for jj in range( cgbead ):
        #        cg_site_of_aa.append(ii)  
    if mode.lower() in ['simplest','compactest','shortest']:
        # expect data to be of the form [(cg_bead_name, n_aa_in_cg, #cg bead)]
        nbead_so_far = 0
        ncg_so_far = 0
        aa_indices_in_cg = []
        cg_site_of_aa = []
        for cgindex,entry in enumerate(shorthand):
            if (isinstance(entry,list) or isinstance(entry,tuple)):
                if len(entry) == 3:
                    cg_bead_name,num_in_cgbead,num_cg_bead = entry
                elif len(entry) == 2:
                    cg_bead_name,num_in_cgbead,num_cg_bead = 'DMY',entry[0],entry[1]
            else:
                raise ValueError('Unrecognized mapping entry for simple shorthand: {}'.format(entry))
            for ii in range(num_cg_bead):
                indices = np.arange(nbead_so_far,nbead_so_far + num_in_cgbead).tolist()
                aa_indices_in_cg.append( (cg_bead_name,indices) ) 
                nbead_so_far += num_in_cgbead
                for jj in range(num_in_cgbead):
                    cg_site_of_aa.append( (cg_bead_name,ncg_so_far) )
                ncg_so_far += 1

    return aa_indices_in_cg, cg_site_of_aa


def generate_single_mapping(mode,extra=None):
    '''
    Generate mapping for a single chain

    Parameters
    ----------
    mode : str or None
        can be a [#atoms x 1] array to map from atom to CG site (i.e. cg site of atom i). Each entry should be tuple (cgbeadname,identifier)
    extra : several
        a mapping array of one of the accepted formats, or a topology for use with 121 or 1res21

    Return
    ------
    aa_indices_in_cg : list of lists
        aa_indices_in_cg[i] returns list of aa indices in cg bead i
    chain_indices_of_cg : list
        chain_id[i] returns the chain index of CG bead i. Hopefully is contiguous...
    '''
    if type(mode) is str:
        if isinstance(extra,mdtraj.Topology):
            top = extra
            if mode.lower() in ['121','1to1']:
                print('doing one-to-one mapping of AA bead to CG bead')
                cg_bead_names_len = top.n_atoms
                aa_indices_in_cg = [ (a.name,[a.index]) for a in top.atoms ]
                cg_site_of_aa = [ (a.name,a.index) for a in top.atoms ]
            elif mode.lower() in ['1res21','1resto1']:
                print('doing one residue to one CG bead')
                aa_indices_in_cg = [ (r.name,[a.index for a in r.atoms]) for r in top.residues ]
                cg_site_of_aa = [ (a.residue.name,a.residue.index) for a in top.atoms ]
        elif mode.lower() in ['cg_site_of_aa']:
            #if len(list(extra)) != top.n_atoms:
            #    raise ValueError('Fed in mapping array of CG site for each atom, should have dimensions #atoms x 1')
            print('using full custom mapping, extra input is [#atoms x 1] array specifying the CG site ID of each atom')
            cg_site_of_aa = extra
            cgsites = [] #ordered in sequence of occurrence
            aa_indices_in_cg = []
            for cgsite in extra: #this is only really efficient if a single molecule, since it doesn't assume contiguity.
                cgsite_identifier = cgsite[1]
                cgsite_name = cgsite[0]
                if cgsite_identifier not in cgsites:
                    cgsites.append(cgsite_identifier)
                    indices_in_cg = [ ii for ii,siteID in enumerate(extra) if siteID[1] == cgsite_identifier ]
                    aa_indices_in_cg.append( (cgsite_name,indices_in_cg) )
        elif mode.lower() in ['aa_indices_in_cg']:
            aa_indices_in_cg = extra
            num_aa = np.array([ len( entry[1] ) for entry in aa_indices_in_cg ]).sum()

            tmp = [[]] * num_aa
            aa_indices = np.zeros( num_aa )
            nbead_so_far = 0
            for isite, site_entry in enumerate(aa_indices_in_cg):
                for index in site_entry[1]:
                    tmp[nbead_so_far] = (site_entry[0],isite)
                    aa_indices[nbead_so_far] = index
                    nbead_so_far += 1
            #next, need to sort
            sorted_indices = np.argsort( aa_indices )
            cg_site_of_aa = [ tmp[index] for index in sorted_indices ]
             

        elif mode.lower() in ['simple','contiguous']:
            print('using compact mapping mode assuming contiguous sequencing of AA beads. Extra input should be array of # AA beads in each CG bead')
            aa_indices_in_cg, cg_site_of_aa = generate_single_mapping_shorthand('simple',extra)
            '''
            aa_indices_in_cg = []
            chain_index_of_cg = []
            nbead_so_far = 0
            for num_in_cgbead in extra:
                indices_in_cg = list(np.arange(nbead_so_far,nbead_so_far + num_in_cgbead)) 
                aa_indices_in_cg.append( indices_in_cg )
                nbead_so_far += num_in_cgbead
            '''

        elif mode.lower() in ['simplest','compactest']:
            print('using most compact mapping mode assuming contiguous sequencing of AA beads. Extra input should be array of # AA beads in each CG bead')
            aa_indices_in_cg, cg_site_of_aa = generate_single_mapping_shorthand('simplest',extra)

        #elif isinstance(mode,(np.ndarray,list)):
    else:
        raise ValueError('non-string mode not implemented yet')

    #sanity check. This function meant to be used on a single chain (molecule), so all the chain indices better be the same one!
    if isinstance(extra,mdtraj.Topology):
        chain_indices_of_cg = []
        for cg_index,entry in enumerate(aa_indices_in_cg):
            indices_in_cg = entry[1]
            chain_indices_in_cgbead = [ top.atom(ii).residue.chain.index for ii in indices_in_cg ] 
            if all( chain_index == chain_indices_in_cgbead[0] for chain_index in chain_indices_in_cgbead ):
                chain_indices_of_cg.append( chain_indices_in_cgbead[0] )
            else:
                 raise ValueError('chain indices {} of atoms in cg bead {} span more than one chain!'.format(chain_indices_in_cgbead,cg_index))
        
    #return aa_indices_in_cg, residue_indices_of_cg, chain_indices_of_cg
    return aa_indices_in_cg, cg_site_of_aa


def generate_system_mapping(mappings):
    ''' 
    Generate the AA->CG mapping for an entire system, from single-chain mappings

    Parameters
    ----------
    mapping: list of tuples
        each tuple should be (aa_index_in_cg_mapping_of_molecule, # occurrences of molecule)

    Return
    ------
    aa_indices_in_cg : list of lists
        aa_indices_in_cg[i] returns list of aa indices in cg bead i

    Notes:
    ------
    The sequence of mapping should correspond to that of the AA system!
    Assume each *molecule* is its own new chain.
    In the future, would be nice to have it s.t. just need to define mapping at residue level, knows how to handle patches, and can get the chain-level mapping. This will necessitate knowing how to manipulate the xml topology definitions.
    '''
    aa_indices_in_cg = []
    n_atoms_total = 0
    n_chains_total = 0
    
    for mapping in mappings:

        molecule_mapping = mapping[0]
        n_replicates = mapping[1] #of replicates
        n_beads_in_mapping = np.sum([ len(cgbead_mapping) for cgbead_mapping in molecule_mapping ])

        cg_bead_names = [cgbead_mapping[0] for cgbead_mapping in molecule_mapping]
        for ii in range(n_replicates):
            new_indices = [ list( np.array(cgbead_mapping[1]) + n_atoms_total ) for cgbead_mapping in molecule_mapping ]

            aa_indices_in_cg.append( list(zip(cg_bead_names,new_indices)) )  #this way, can later iterate through by molecule/chain!
            n_atoms_total += n_beads_in_mapping

    return system_aa_indices_in_cg


#def map(top,mapping_AA2CG,beadtypes):
    '''
    Parameters
    ----------
    top : mdtraj topology
    mapping : array
        following the output of my generate() helper functions, i.e. returns [ [aa_indices_in_cgbead1], [aa_indices_in_cgbead2], ... ], grouped by chain

    Return
    ------
    newtop : mdtraj topology

    Notes
    -----
    Current implementation best for mapping single chain/small systems. For larger systems, will have to do a parallel implementation

    Nice thing is that this code will figure out the bonds automatically!
    Only check for bonding if intrachain.
    '''
#    for ic in chains:
#        for ir in residues: #I don't use residues in CG model?
#            for ia in atomlistofresidue:
#                add new atom
#

# ===== Processing, Mapping and Synthesizing =====
def process_pdbfile(filename):
    '''
    For if the mapping uses the residue name column of pdbfile as the cg bead
    '''
    t = mdtraj.load(filename)
    cg_site_of_aa = [(a.residue.name,a.residue.index) for a in t.top.atoms]
    aa_indices_in_cg,_ = generate_single_mapping(mode='cg_site_of_aa',extra=cg_site_of_aa)

    return aa_indices_in_cg, cg_site_of_aa

def process_mappingfile(filename,mode,customfield=None):
    '''
    For if the (single chain) mapping uses a yaml/text-based specification of the array
    '''
    with open(filename,'r') as stream:
        mapping = yaml.load(stream)
   
    if mode.lower() == 'pdb':
        my_mapping = mapping['pdbfile_mapping']
        aa_indices_in_cg, cg_site_of_aa = process_pdbfile( my_mapping )
        shorthand=None
    elif mode.lower() in ['aa','aa_indices_in_cg','aa_index','aa_indices']:
        my_mapping = mapping['aa_indices_in_cg']
        aa_indices_in_cg, cg_site_of_aa = generate_single_mapping( mode='aa_indices_in_cg',extra=my_mapping )
        shorthand=None
    elif mode.lower() in ['cg','cg_site_of_aa','cg_site','cg_sites']:
        my_mapping = mapping['cg_site_of_aa']
        aa_indices_in_cg, cg_site_of_aa = generate_single_mapping( mode='cg_site_of_aa',extra=my_mapping )
        shorthand=None
    elif mode.lower() in ['short','shorthand','simple','contiguous']:
        if customfield is None:
            my_mapping = mapping['shorthand']
        else:
            my_mapping = mapping[customfield]
        aa_indices_in_cg, cg_site_of_aa = generate_single_mapping_shorthand( mode='short',shorthand=my_mapping )
        shorthand = my_mapping
    elif mode.lower() in ['shortest','shortesthand','simplest']:
        if customfield is None:
            my_mapping = mapping['shorthand']
        else:
            my_mapping = mapping[customfield]
        aa_indices_in_cg, cg_site_of_aa = generate_single_mapping_shorthand( mode='shortest',shorthand=my_mapping)
        shorthand = my_mapping
    else:
        print('caution, trying to use a custom mode {} that may not be implemented yet').format(mode)
        if customfield is None:
            my_mapping = mapping[mode]
        else:
            my_mapping = mapping[customfield]
        aa_indices_in_cg, cg_site_of_aa = generate_single_mapping_shorthand( mode=mode,extra=my_mapping)
        shorthand = None
        #raise ValueError('specified mode has not been implemented yet')
       
    return aa_indices_in_cg, cg_site_of_aa, shorthand

def map_single(traj,mapping):
    '''
    take in a mapping array, and generate a new mapped topology for the chain, retaining bonds!
    Parameters
    ----------
    traj : mdtraj trajectory
    mapping : array
        should be aa_indices_in_cgbead array format, of a single chain
    '''
    cg_top = mdtraj.Topology()
    #=== first add all the species ===
    print('creating topology: adding species')
    ch = cg_top.add_chain()
    n_beads = len(mapping)
    xyz = np.zeros([n_beads,3])
    cg_beads = []
    for ic,cgbead_entry in enumerate(mapping):
        cgbead_name,aa_indices = cgbead_entry
        resname = traj.top.atom(aa_indices[0]).residue.name #use first atom's residue as the resname. alternative: use most common occurrence
        masses = np.array([traj.top.atom(ia).element.mass for ia in aa_indices])
        com = np.sum(traj.xyz[0,aa_indices,:] * masses[:,None],0) / masses.sum()
        #print( mdtraj.compute_center_of_mass(traj.atom_slice(aa_indices))[0] )
        #print( com )
        xyz[ic,:] = com

        res = cg_top.add_residue( resname,ch )
        element = mdtraj.element.virtual #for now... later on can create elements for CG simulation from cg bead name
        bead = cg_top.add_atom(cgbead_name, element, res)
        cg_beads.append(bead)

    #=== figure out connectivity ===
    #loop over each cg bead, then over other cg beads, then check bonding
    print('creating topology: figuring out bonding')
    bondgraph = traj.top.to_bondgraph()
    for ii in range(n_beads):
        aa_indices1 = mapping[ii][1]
        for jj in range(ii+1,n_beads):
            aa_indices2 = mapping[jj][1]
            bonded = any([ bondgraph.has_edge(traj.top.atom(a1),traj.top.atom(a2)) for a1 in aa_indices1 for a2 in aa_indices2 ])
            if bonded:
                cg_top.add_bond(cg_beads[ii],cg_beads[jj])

    #=== create mapped trajectory and save ===
    cg_traj = mdtraj.Trajectory([xyz],cg_top)

    return cg_traj




# ===== Tests =====
#t = mdtraj.load('SDS.pdb')
#try to figure out the mapping
#rando = [(a.residue.name,a.residue.index) for a in t.top.atoms]
#print('start test')
#print(rando)
#chain_mapping,cg_indices = generate_single_mapping(t.topology,mode='full',extra=rando)

#t = mdtraj.load('pS20.pdb')
#rando = np.random.randint(0,5,160)
#generate_system_mapping( [(chain_mapping,5)] )
#chain_mapping2,cg_indices2 = generate_single_mapping(t.topology,mode='simple',extra=[8]*t.n_residues)


# Experiment 2: map from an annotated pdb file
'''
filename = 'SDS_ordered_2gro.pdb'
t = mdtraj.load(filename)
prefix = os.path.splitext(filename)[0]
aa_indices_in_cg, cg_site_of_aa = process_pdbfile('SDS.pdb')
traj_mapped = map_single(t,aa_indices_in_cg)
traj_mapping = generate_pdb_mapping(t, cg_site_of_aa)
save(prefix, aa_indices_in_cg, cg_site_of_aa, traj_mapped = traj_mapped, traj_mapping = traj_mapping )
'''


# Experiment 3: map from a given mapping array
filename = 'SDS_ordered_2gro.pdb'
filemap = 'SDS.yaml'
t = mdtraj.load(filename)

## Test pdbfile
print('--- test0 ---')
prefix,mode = 't0','pdb'
aa_indices_in_cg, cg_site_of_aa, _ = process_mappingfile(filemap,mode)
traj_mapped = map_single(t,aa_indices_in_cg)
traj_mapping = generate_pdb_mapping(t, cg_site_of_aa)
save(prefix, aa_indices_in_cg, cg_site_of_aa, traj_mapped = traj_mapped, traj_mapping = traj_mapping )

## Test aa_indices_in_cg
print('--- test1 ---')
prefix,mode = 't1','aa'
aa_indices_in_cg, cg_site_of_aa, _ = process_mappingfile(filemap,mode)
traj_mapped = map_single(t,aa_indices_in_cg)
traj_mapping = generate_pdb_mapping(t, cg_site_of_aa)
save(prefix, aa_indices_in_cg, cg_site_of_aa, traj_mapped = traj_mapped, traj_mapping = traj_mapping )

## Test cg_site_of_aa
print('--- test2 ---')
prefix,mode = 't2','cg'
aa_indices_in_cg, cg_site_of_aa, _ = process_mappingfile(filemap,mode)
traj_mapped = map_single(t,aa_indices_in_cg)
traj_mapping = generate_pdb_mapping(t, cg_site_of_aa)
save(prefix, aa_indices_in_cg, cg_site_of_aa, traj_mapped = traj_mapped, traj_mapping = traj_mapping )

## Test shorthand
print('--- test3 ---')
prefix,mode = 't3','short'
aa_indices_in_cg, cg_site_of_aa, shorthand = process_mappingfile(filemap,mode)
traj_mapped = map_single(t,aa_indices_in_cg)
traj_mapping = generate_pdb_mapping(t, cg_site_of_aa)
save(prefix, aa_indices_in_cg, cg_site_of_aa, traj_mapped = traj_mapped, traj_mapping = traj_mapping, shorthand=shorthand )

## Test shortesthand
print('--- test4 ---')
prefix,mode = 't4','shortest'
#aa_indices_in_cg, cg_site_of_aa = process_mappingfile(filemap,mode)
aa_indices_in_cg, cg_site_of_aa, shorthand = process_mappingfile(filemap,mode,customfield='shortesthand')
traj_mapped = map_single(t,aa_indices_in_cg)
traj_mapping = generate_pdb_mapping(t, cg_site_of_aa)
save(prefix, aa_indices_in_cg, cg_site_of_aa, traj_mapped = traj_mapped, traj_mapping = traj_mapping, shorthand=shorthand)
