#!/usr/bin/env python3
"""
mican-analysis skill — environment setup script.

Run this script once before using the skill:
    python setup_env.py

Checks for pymican and installs it if missing.
pymican requires a C compiler (gcc) to build the MICAN binary.
"""

import sys
import subprocess
import shutil


def check_pymican() -> bool:
    """Return True if pymican is importable and the binary works."""
    try:
        from pymican import mican, BINFILEPATH
        import pathlib
        return pathlib.Path(BINFILEPATH).exists()
    except Exception:
        return False


def check_compiler() -> bool:
    """Return True if gcc or cc is available."""
    return shutil.which('gcc') is not None or shutil.which('cc') is not None


def install_pymican():
    """Install pymican via pip."""
    print("Installing pymican...")
    result = subprocess.run(
        [sys.executable, '-m', 'pip', 'install', 'pymican'],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print("pip install failed:")
        print(result.stderr[-500:])
        return False
    print("pymican installed successfully.")
    return True


def verify():
    """Verify that pymican works end-to-end."""
    try:
        from pymican import mican
        m = mican()
        print("pymican import: OK")
        print(f"MICAN binary:   OK")
        return True
    except Exception as e:
        print(f"Verification failed: {e}")
        return False


def main():
    print("=" * 50)
    print("mican-analysis skill: environment setup")
    print("=" * 50)

    # 1. Already installed?
    if check_pymican():
        print("✅ pymican is already installed and working.")
        verify()
        return

    # 2. Compiler check
    if not check_compiler():
        print("⚠ No C compiler found (gcc/cc).")
        print("  pymican requires gcc to compile the MICAN binary.")
        print("  Install gcc first:")
        print("    macOS:  xcode-select --install  (or brew install gcc)")
        print("    Ubuntu: sudo apt install gcc")
        print("    CentOS: sudo yum install gcc")
        sys.exit(1)

    # 3. Install
    if not install_pymican():
        print("\n⚠ Installation failed.")
        print("  You can also build from source:")
        print("    git clone https://github.com/ShintaroMinami/pymican")
        print("    cd pymican && make && python setup.py install")
        sys.exit(1)

    # 4. Verify
    if verify():
        print("\n✅ Setup complete. The mican-analysis skill is ready to use.")
    else:
        print("\n⚠ Installation succeeded but verification failed.")
        print("  Try restarting your Python session and running the skill again.")
        sys.exit(1)


if __name__ == '__main__':
    main()
