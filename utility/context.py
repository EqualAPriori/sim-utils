import sys,os,subprocess
sys.path.insert(1, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import sim


def get_path():
  return 'python path:\n    {}'.format(sys.path)

def getGitBranch(path):                                                
    """ Checks whether working directory (of .blend file) contains a git repository.
        Returns branch if repository is found.                           
        
        Notes:
          https://www.programcreek.com/python/?CodeExample=get+git+branch  
    """
    try:
      output = str(subprocess.check_output(['git', 'branch'],        
                   cwd=path, universal_newlines=True)
                  )
      branch = [a for a in output.split('\n') if a.find('*') >= 0][0]
      return branch[branch.find('*') + 2 :] 
    except subprocess.CalledProcessError:
      print('ProcessError') 
      return None 
    except FileNotFoundError:
      log("No git repository found.", "ERROR") 
    return None                                                          

def get_devinfo(module):
  ''' To document the version of a module we're running with

  Args:
    module (module): the module to inspect/document

  Returns:
    output (str)

  Todo:
    once I have version numbers, can use git describe to get slightly more informative information
  '''
  path = os.path.dirname(module.__file__)
  whitespace = ' ' * len(module.__name__)
  s = '... {} repo:           {}'.format(module.__name__, path)

  # --- branch name with `git rev-parse --abbrev-ref HEAD`
  #print('... {} branch:\t{}'.format(whitespace, getGitBranch(simpath) ))            
  branchname = subprocess.check_output(['git','rev-parse','--abbrev-ref','HEAD'], cwd=path)
  s = s + '\n' + '... {} branch:         {}'.format(whitespace,branchname)

  # --- get commit info with git show -s --date=short --pretty='format:%h (%s, %ad)'
  #print('... {} commit:\t{}'.format(whitespace, subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=simpath)))
  #commitinfo = subprocess.check_output("git show -s --date=short --pretty='format:%h (%s, %ad)'".split(), cwd=path)
  #git show -s --pretty=format:"commit: %h%nauthor: %an%ndate:   %ad%nmsg:    %<(80,trunc)%s"
  #commitinfo = subprocess.check_output(["git","show","-s",'--pretty=format:hash:   %h%nauthor: %an%ndate:   %ad%nmsg:    %<(80,trunc)%s'],cwd=path)
  commitinfo = subprocess.check_output(["git","show","-s",'--pretty=format:%h    %an    %ad%'],cwd=path)
  commitmsg  = subprocess.check_output(["git","show","-s",'--pretty=format:%<(80,trunc)%s'],cwd=path)
  s = s + '... {} commit:         {}'.format(whitespace, commitinfo)
  s = s + '\n' + '... {} message:        {}'.format(whitespace, commitmsg)

  # --- get uncommitted changes to watch out for
  uncommited_files = subprocess.check_output(['git', 'diff', '--name-only', 'HEAD'], cwd=path)
  s = s + '\n' + "... {} uncommitted:\n{}".format(whitespace, uncommited_files)

  # --- optional
  #print('... description:\t{}'.format(subprocess.check_output(['git', 'describe','--always'])))
  #print('... cwd:\t{}'.format(os.getcwd()))                              

  #print(s)
  return(s)
      




