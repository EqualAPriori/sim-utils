### Helper classes and functions to create force fields
import sim
from collections import namedtuple
import system
from dotdict import dotdict

### ===== Wishlist =====
#   Bonding
#   Spline

### ===== DataTypes (for holding parameters), defined to enable dot notation) =====
#   alternative is to use dictionary and key notation
class Parameter:
    def __init__(self,value,fixed=False):
        self.val = value
        self.fixed = fixed
    def __repr__(self):
        return '{{value: {}, fixed: {}}}'.format(self.val,self.fixed)
    #def __str__(self):
    #    return 'value: {}, fixed: {}'.format(self.val,self.fixed)


#Attempt at creating data class that allows the: obj.param and obj.param.fixed notation... actually not sure how to get obj.param to return the value, but obj.param.fixed to return the boolean value for whether or not the parameter is fixed...
#class LJG:
#    def __init__(self,B=(1.0,False),Kappa=(0.25,True),Dist0=(0.0,True),Sigma=(1.0,True),Epsilon=(0.0,True),Label=None):
#
#ah... namedtuple w/defaults is only for python 3.6+.
#
#LJG = namedtuple('LJG','B Kappa Dist0 Sigma Epsilon Label', 
#        defaults=[Parameter(1.0,False),Parameter(0.25,True),Parameter(0.0,True),Parameter(1.0,True),Parameter(0.0,True),None])
#ExtSin = namedtuple('ExtSin','UConst NPeriods PlaneAxis PlaneLoc Label',
#        defaults=[Parameter(0.0,True),Parameter(1.,True),Parameter(0,True),Parameter(0.,True),None])
#
# Ulimately decide the parameter library's structure:
# - type (e.g. ljg, ewald, etc)
# - label
# - cutoff
# - other pertinent parameters
# - key-value pairs for all the other parameters, stored as Parameter objects
#
# alternative implementation: split off values and fixed into two separate dictionaries/arrays
# - values = dict
# - fixed = dict
#
#and I want the input for params and fixed to be... dicts?
#
def set_parameters(targetdict, params):
    """Helper function for the more specific functions for setting different potential types' parameters
    Notes
    -----
    Allows for 4 different notations:
    * bool only: interpreted specifying fixed/not
    * int/float: interpreted as value
    * list/tuple: interpreted as (value,fixed) pair
    # dict: interpreted as kwarg pairs for Parameter
    """
    if isinstance(params,dict):
        for key,value in params.items():
            if isinstance(value,(bool)):
                targetdict[key].fixed = value
            elif isinstance(value,(int,float)):
                targetdict[key].val = value
            elif isinstance(value,(list,tuple)):
                targetdict[key] = Parameter(*value)
            elif isinstance(value,dict):
                targetdict[key] = Parameter(**value)
            else:
                raise ValueError('Unsuported value type for parameter {}'.format(key))
    else:
        raise ValueError('params must be input as a dictionary')

def param_ljg(params, cut, fixed_default=True, label=None):
    """Set up a LJG dictionary
    Notes
    -----
    Caution, no type checking implemented
    """
    tmp = {'type':'ljg','Label':label,'Cut':cut,'B':Parameter(1.0,fixed_default), 'Kappa':Parameter(0.25,fixed_default), 'Dist0':Parameter(0.0,fixed_default), 'Sigma':Parameter(1.0,fixed_default), 'Epsilon':Parameter(0.0,fixed_default)}
    set_parameters(tmp,params)
    """
    if isinstance(fixed,dict):
        for key,value in fixed.items():
            tmp[key].fixed = value
    elif fixed is None:
        pass 
    else:
        raise ValueError('fixed parameter must be a dictionary')
    """
    return dotdict(tmp)

def param_ext_sin(params, fixed_default=True, label=None):
    tmp = {'type':'ext_sin','Label':label, 'UConst':Parameter(0.0,True),'NPeriods':Parameter(1.,True),'PlaneAxis':Parameter(0,True),'PlaneLoc':Parameter(0.,True)}
    set_parameters(tmp,params)
    return dotdict(tmp)

def param_sm_coul_ew_corr(params, cut, sm_coul_shift=True, fixed_default=True, label=None):
    tmp = {'type':'sm_coul_ew_corr','Label':label,'Cut':cut,'Shift':sm_coul_shift,'BornA':Parameter(1.0,fixed_default),'Coef':Parameter(1.0,fixed_default)}
    set_parameters(tmp,params)
    return dotdict(tmp)

def param_ewald(params, cut, exclude_bond_ord=0, shift=True, label=None):
    tmp = {'type':'ewald','Label':label, 'Cut':cut, 'Shift':shift, 'ExcludeBondOrd':exclude_bond_ord, 'Coef':Parameter(0.0,True)}
    set_parameters(tmp,params)
    return dotdict(tmp)

def param_bond(params, fixed_default=True, label=None):
    tmp = {'type':'bond','Label':label, 'Dist0':Parameter(0.0,fixed_default),'FConst':Parameter(1.0,fixed_default)}
    set_parameters(tmp,params)
    return dotdict(tmp)


### ===== Set up the force field =====
#   Right now assumes fairly simple tuples/indices, i.e. does not yet interpret more complicated pair filters of the form ( (atoms), (atoms) )

def create_forcefield(Sys, top, interactions, special_options=None):
    sim_atom_types = { atom.Name:atom for atom in Sys.World.AtomTypes }
    filters={} 
    force_field = []
    if special_options is None:
        special_options = system.make_special_options()
    def atomtypes_in_system(atype_names):
        return all( [atype_name in sim_atom_types for atype_name in atype_names] )

    for atom_tuple,interaction_list in interactions.items():
        for params in interaction_list:
            #... Lennard Jones Gaussian ...
            if params['type'].lower() in ['ljg']:
                if atomtypes_in_system(atom_tuple):
                    f = sim.atomselect.PolyFilter( [sim_atom_types[atom_tuple[0]],sim_atom_types[atom_tuple[1]]] )
                    filters[atom_tuple] = f
                    print('Adding LJG for {}'.format(atom_tuple))
                    p = sim.potential.LJGaussian( Sys, Filter=f, Cut=params['Cut'], Fixed=True,
                            B=params.B.val, Kappa=params.Kappa.val, Dist0=params.Dist0.val, Sigma=params.Sigma.val, Epsilon=params.Epsilon.val, Label=params.Label )

                    p.Param.B.Fixed = params.B.fixed
                    p.Param.B.Min = -100. #default, don't restrict positive
                    p.Param.Kappa.Fixed = params.Kappa.fixed
                    p.Param.Dist0.Fixed = params.Dist0.fixed
                    p.Param.Sigma.Fixed = params.Sigma.fixed
                    p.Param.Epsilon.Fixed = params.Epsilon.fixed
                    force_field.append(p)
                else:
                    print('ATOM TYPES in {} NOT PRESENT, SKIPPING').format(atom_tuple)

            #... Smeared Coulomb ...
            if params['type'].lower() in ['sm_coul_ew_corr']:
                if special_options['neutralize']:
                    print('Forcing system to be neutral, skipping ewald interaction')
                else:
                    if atomtypes_in_system(atom_tuple):
                        f = sim.atomselect.PolyFilter( [sim_atom_types[atom_tuple[0]],sim_atom_types[atom_tuple[1]]] )
                        filters[atom_tuple] = f
                        print('Adding smeared coulomb ewald correction for {}'.format(atom_tuple))
                        p = sim.potential.SmearedCoulombEwCorr( Sys, Filter=f, Cut=params.Cut, Shift=params.Shift, Coef=params.Coef.val, FixedCoef=params.Coef.fixed, BornA=params.BornA.val, FixedBornA=params.BornA.fixed, Label=params.Label )
                        force_field.append(p)
                    else:
                        print('ATOM TYPES in {} NOT PRESENT, SKIPPING').format(atom_tuple)

            #... Ewald ...
            if params['type'].lower() in ['ewald']:
                if special_options['neutralize']:
                    print('Forcing system to be neutral, skipping ewald interaction')
                else:
                    print('Adding ewald correction for {}'.format(atom_tuple))
                    p = sim.potential.Ewald( Sys, ExcludeBondOrd=params.ExcludeBondOrd, Cut=params.Cut, Shift=params.Shift, Coef=params.Coef.val, FixedCoef=params.Coef.fixed, Label=params.Label )
                    force_field.append(p)


            #... External Sine ...
            if params['type'].lower() in ['extsin','ext_sin','external_sinusoid','external_sine']: #External Sinusoid
                #interpret tuple as a single filter with many atoms, rather than a pair filter
                #f = sim.atomselect.PolyFilter( sim.atomselect.Filter([sim_atom_types[atom_name] for atom_name in atom_tuple]) )
                if atomtypes_in_system(atom_tuple):
                    f = sim.atomselect.Filter([ sim_atom_types[atom_name] for atom_name in atom_tuple ]) 
                    filters[atom_tuple] = f
                    print('Adding external sinusoid for {}'.format(atom_tuple))
                    p = sim.potential.ExternalSinusoid(Sys, Filter=f, Fixed=True, UConst=params.UConst.val, NPeriods=params.NPeriods.val, PlaneAxis=params.PlaneAxis.val, PlaneLoc=params.PlaneLoc.val, Label=params.Label)
                    p.Param.UConst.Fixed = params.UConst.fixed
                    p.Param.NPeriods.Fixed = params.NPeriods.fixed
                    force_field.append(p)
                else:
                    print('ATOM TYPES in {} NOT PRESENT, SKIPPING').format(atom_tuple)


            #... Harmonic Bonds ...
            if params['type'].lower() in ['bond','harmonic','harmonic_bond','harmonicbond']:
                if atomtypes_in_system(atom_tuple):
                    print('Adding harmonic bond for {}'.format(atom_tuple))
                    f = sim.atomselect.PolyFilter( [sim_atom_types[atom_tuple[0]],sim_atom_types[atom_tuple[1]]], Bonded=True )
                    p = sim.potential.Bond(Sys, Filter=f, Fixed=True, Dist0=params.Dist0.val, FConst=params.FConst.val, Label=params.Label)
                    p.Param.Dist0.Fixed = params.Dist0.fixed
                    p.Param.FConst.Fixed = params.FConst.fixed
                    force_field.append(p)
                else:
                    print('ATOM TYPES in {} NOT PRESENT, SKIPPING').format(atom_tuple)

            print('    {}'.format(params))
    return force_field


def set_param_from_file(force_field, ff_file):
    with open(ff_file, 'r') as of: s = of.read()
    print('... Setting FF with file {} ...\n{}'.format(ff_file,s))
    force_field.SetParamString(s)      


def update_specific(force_field, params):
    """Update specific terms of an existing Sys.ForceField
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


def get_free_parameters(force_field):
    free_params = []
    for potential in force_field:
        for ip,param in enumerate(potential.Param.Names):
            fixed = potential.Fixed[ip]
            if not fixed:
                free_params.append( (potential.Name,param) )
    return free_params


def toggle_fixed(force_field,params,fixed):
    ffdict = {p.Name:p for p in force_field}
    for param in params:
        potential_name = param[0]
        parameter_name = param[1]
        ffdict[potential_name].__getattr__(parameter_name).Fixed = fixed



