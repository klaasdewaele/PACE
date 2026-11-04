"""
setup.py — filesystem setup utilities for the PACE pipeline.

Provides helper functions for validating the input data directory and
creating a uniquely named session output folder for each pipeline run.
"""

import logging
import os
from datetime import datetime

def input_setup(input_path):
    '''Check if input_path exists. If not, exit.'''
    if os.path.exists(input_path):
        logging.info(f'Input directory: {input_path}')
        return input_path
    else:
        logging.warning(f'No valid input directory provided. Quitting...')
        quit()

def create_session_folder(session_folder):
    '''Recursive function to avoid duplicate session_folder when executing script parallel at same time.'''
    try:
        os.mkdir(session_folder)
        logging.info(f'Created output folder at {session_folder}')
        return session_folder
    except FileExistsError:
        logging.debug(f'FileExistsError {session_folder}: adapting path name')
        seconds_plus_one = str(int(session_folder.split('-')[-1]) + 1)
        new_session_folder = '-'.join(session_folder.split('-')[:-1] + [seconds_plus_one])
        session_folder = create_session_folder(new_session_folder)
        return session_folder

def output_setup(output_path, home, label):
    '''Creates output directory; if not provided by user, an output directory is created in a home directory defined in main.py.'''
    session_time = datetime.now().strftime('%d-%m-%Y-%H-%M-%S')
    if os.path.exists(output_path):
        session_folder_final = os.path.join(output_path, f'{label}_{session_time}')
        os.mkdir(session_folder_final)
        logging.info(f'Session folder created in provided output folder: {output_path}')
    else:
        logging.info('No valid output directory provided. Output folder created in home directory.')
        session_folder = os.path.join(home, f'{label}_{session_time}')
        session_folder_final = create_session_folder(session_folder)
    return session_folder_final, session_time
