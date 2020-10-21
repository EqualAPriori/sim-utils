# Minimal dotdict implementation, not fully featured
#   most critically, get() returns None if key doesn't exist instead of throwing an attribute error, be warned!!!
#   adapted from: https://stackoverflow.com/questions/2352181/how-to-use-a-dot-to-access-members-of-dictionary/28463329
# for full implementation, see DotMap (https://github.com/drgrib/dotmap)
# also described on above stackoverflow link
class dotdict(dict):
    """
    A dictionary supporting dot notation.
    """
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __init__(self, *args, **kwargs):
        #print(args)
        #print(kwargs)
        #super(dict,self).__init__(*args, **kwargs) #python 3+ don't need explicit super() arguments
        self.update(*args,**kwargs)
        #print(self)
        for k, v in self.iteritems(): #python 3+ use items()
            if isinstance(v, dict):
                self[k] = dotdict(v)

    def lookup(self, dotkey):
        """
        Lookup value in a nested structure with a single key, e.g. "a.b.c"
        """
        path = list(reversed(dotkey.split(".")))
        v = self
        while path:
            key = path.pop()
            if isinstance(v, dict):
                v = v[key]
            elif isinstance(v, list):
                v = v[int(key)]
            else:
                raise KeyError(key)
        return v

if __name__ == '__main__':
     print('... Testing dotdict ...')
     dd = dotdict({"a": 1, "b": {"c": "# it works", "d": [1, 2, {"e": 123}]}})
     print(dd.b.c)
     print("#", isinstance(dd, dict))
     print("#", dd.lookup("b.d"))
     # it works
     # True
     # [1, 2, {'e': 123}]
