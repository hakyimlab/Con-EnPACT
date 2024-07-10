import os,sys
import subprocess

def make_enpact_prediction(input_epigenome, output_predictions, model_path, script_directory):

    com = [
        "Rscript",
        "--vanilla",
        os.path.join(script_directory,"predict_from_epigenome.R"),
        f"--input_file={input_epigenome}",
        f"--output_file={output_predictions}",
        f"--trained_model={model_path}"
    ]

    print(" ".join(com))

    subprocess.run(com)