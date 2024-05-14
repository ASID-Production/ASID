from Source import MainWindow
import sys
import os
import logging
from datetime import datetime

time = datetime.today().strftime('%Y-%m-%d')

dirs = ['temp', 'logs']
for dir in dirs:
    if not os.path.exists(f'{dir}/'):
        os.mkdir(dir)

if '--debug' in sys.argv or '-d' in sys.argv:
    path = os.path.join(os.path.dirname(__file__), f'logs/{time}-debug.log')
    logging.basicConfig(level=logging.INFO, filename=path, filemode='a', format="%(asctime)s - [%(levelname)s] -  %(name)s - (%(filename)s).%(funcName)s(%(lineno)d) - %(message)s")
else:
    path = os.path.join(os.path.dirname(__file__), f'logs/{time}-error.log')
    logging.basicConfig(level=logging.WARNING, filename=path, filemode='a', format="%(asctime)s - [%(levelname)s] -  %(name)s - (%(filename)s).%(funcName)s(%(lineno)d) - %(message)s")
module_dir = os.path.normpath(os.path.join(os.path.split(__file__)[0], '../module'))
base_dir = os.path.normpath(os.path.split(__file__)[0])
sys.path.append(module_dir)
sys.path.append(base_dir)

try:
    MainWindow.show()
except Exception as ex:
    logger = logging.getLogger(__name__)
    logger.exception('Exception')
    raise ex
