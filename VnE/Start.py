import sys
import os
base_dir = os.path.normpath(os.path.split(__file__)[0])
sys.path.append(base_dir)

from Source import MainWindow
import logging
from datetime import datetime
import ctypes
import os
if os.name == 'nt':
    myappid = u'ASID-team.ASID.ASID-dev.1' # arbitrary string
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

time = datetime.today().strftime('%Y-%m-%d')
dirs = ['temp', 'logs']
for dir in dirs:
    if not os.path.exists(f'{os.path.dirname(__file__)}/{dir}/'):
        os.mkdir(f'{os.path.dirname(__file__)}/{dir}')

if '--debug' in sys.argv or '-d' in sys.argv:
    path = os.path.normpath(os.path.join(os.path.dirname(__file__), f'{os.path.dirname(__file__)}/logs/{time}-debug.log'))
    logging.basicConfig(level=logging.INFO, filename=path, filemode='a', format="%(asctime)s - [%(levelname)s] -  %(name)s - (%(filename)s).%(funcName)s(%(lineno)d) - %(message)s")
else:
    path = os.path.normpath(os.path.join(os.path.dirname(__file__), f'{os.path.dirname(__file__)}/logs/{time}-error.log'))
    logging.basicConfig(level=logging.WARNING, filename=path, filemode='a', format="%(asctime)s - [%(levelname)s] -  %(name)s - (%(filename)s).%(funcName)s(%(lineno)d) - %(message)s")
module_dir = os.path.normpath(os.path.join(os.path.split(__file__)[0], '../module'))
sys.path.append(module_dir)

try:
    MainWindow.show()
except Exception as ex:
    logger = logging.getLogger(__name__)
    logger.exception('Exception')
    raise ex
