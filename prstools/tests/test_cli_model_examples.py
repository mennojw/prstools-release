import tempfile
import subprocess
import os, re
import pytest
from prstools.models import PRSCS2 as PRSCS2test # Replace with the actual class

def extract_commands_from_epilog(model_cls):
    lines = model_cls._get_cli_epilog().splitlines()
    lines = lines[1:]  # drop first line
    return lines #[line.split('#')[0].strip() for line in lines if line.strip()]

# @pytest.mark.slow
def test_model_epilog_examples_run():
    commands = extract_commands_from_epilog(PRSCS2test)

    with tempfile.TemporaryDirectory() as tmpdir:
        cwd = os.getcwd()
        os.chdir(tmpdir)
        try:
            cmd = "\n".join(commands)
            ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
            cmd = ansi_escape.sub('', cmd)
            print(f"$ {cmd}\n", flush=True)
            print('--<=>--'*5, flush=True)            
            with subprocess.Popen(
                cmd, shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1  # line-buffered
            ) as proc:
                for line in proc.stdout:
                    print(line, end='', flush=True)
                proc.wait()
                retcode = proc.returncode
#             result = subprocess.run(
#                 cmd, shell=True,
#                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
#                 text=True
#             )
#             print(f"$ {cmd}")
#             print(result.stdout)
            print('NOW TO ASSERT')
            assert result.returncode == 0, f"Command failed: {cmd}\n{result.stdout}"
        finally:
            os.chdir(cwd)
