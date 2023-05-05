from distutils.core import setup
import glob

setup(name='pyARCiS',
      version='1.0',
      description='python ARCiS interface',
      author='Michiel Min',
      author_email='M.Min@sron.nl',
      packages=[''],
      package_data={'': glob.glob('pyARCiS*.so')},
     )

