import json
import os
from collections import OrderedDict
import datetime


# ===== Utility
defaultname = "z.log.json"

datetimeformat = "%Y-%m-%dT%H-%M-%S"
dateformat = "%Y-%m-%d"
def now(obj = False):
  """ return current datetime
  Args:
    obj (bool): True --> return datetime object, False --> return string (default)
  """
  if obj:
    return datetime.datetime.now()
  else:
    return datetime.datetime.now().strftime(datetimeformat)

def today(obj = False):
  """ return today's date
  Args:
    obj (bool): True --> return datetime object, False --> return string (default)
  """
  if obj:
    return datetime.date.today()
  else:
    return datetime.date.today().strftime(dateformat)



# ===== Logging
def log(d,filename=defaultname,msg=None):
    """directly writes to log; no intermediate dictionary object returned"""
    d0 = load(filename)

    changes = []
    for k,v in d.items():
        if k in d0:
            if v != d0[k]:
                d0[k] = v
                changes.append("updated {}".format(k))
        else:
            changes.append("added {}".format(k))
            d0[k] = v

    note = "; ".join(changes)
    if msg is not None:
        note = msg + "; " + note
        
    if note != "":
        print(note)
        timestamp = now()
        d0['history'].append((timestamp,note))
    else:
        print("nothing changed in log")

    with open(filename,'w') as f:
        json.dump( d0, f, indent=2 ) 

            
def load(filename=defaultname):
    if os.path.exists(filename):
        with open(filename,'r') as f:
            d0 = json.load(f,object_pairs_hook = OrderedDict)
    else:
        print('Log not found, starting a new log')
        timestamp = now()
        d0 = OrderedDict([("history",[(timestamp,"created")])])
        with open(filename,'w') as f:
            json.dump( d0, f, indent=2 ) 

    return d0
