#!/usr/bin/env python3
import os
import json
import subprocess
import sys

def get_peacock_command(color_name, color_hex):
    """Generate the VS Code command to set Peacock color"""
    return f'code --command "peacock.changeColorToFavorite" --args "{{\\"colorName\\":\\"{color_name}\\",\\"colorHex\\":\\"{color_hex}\\"}}"'

def detect_environment():
    """Detect current environment and return appropriate color"""
    # Check if in Singularity container
    if os.path.exists('/.singularity.d') or os.environ.get('SINGULARITY_CONTAINER'):
        return ("Singularity Container", "#ff79c6")  # Pink
    
    # Check if on HPC (SSH connection)
    if os.environ.get('SSH_CONNECTION') or os.environ.get('SSH_CLIENT'):
        # Check if in interactive session (common HPC indicators)
        if os.environ.get('SLURM_JOB_ID') or os.environ.get('PBS_JOBID'):
            return ("HPC Interactive", "#ff79c6")  # Magenta
        return ("HPC SSH", "#ff79c6")  # Magenta
    
    # Check Python environment
    if hasattr(sys, 'real_prefix') or (hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix):
        venv_name = os.path.basename(sys.prefix)
        return (f"Python venv: {venv_name}", "#f1fa8c")  # Yellow
    
    if os.environ.get('CONDA_DEFAULT_ENV') and os.environ.get('CONDA_DEFAULT_ENV') != 'base':
        return (f"Conda: {os.environ.get('CONDA_DEFAULT_ENV')}", "#ffb86c")  # Orange
    
    # Default local environment
    return ("Local", "#6272a4")  # Purple (Dracula theme purple)

def main():
    env_name, color = detect_environment()
    print(f"Environment detected: {env_name}")
    print(f"Setting Peacock color to: {color}")
    
    # Set VS Code settings
    settings_path = os.path.join(os.path.dirname(__file__), 'settings.json')
    
    try:
        with open(settings_path, 'r') as f:
            settings = json.load(f)
    except:
        settings = {}
    
    # Update Peacock color
    settings['peacock.color'] = color
    
    with open(settings_path, 'w') as f:
        json.dump(settings, f, indent=2)
    
    print(f"Updated VS Code settings with Peacock color: {color}")

if __name__ == "__main__":
    main()