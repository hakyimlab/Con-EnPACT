import os,sys
import subprocess


def make_enpact_predictions(input_file, output_file, model_path, scripts_directory):
    if os.path.exists(output_file):
        print(output_file," already exists")
        return

    com = [
        "Rscript",
        "--vanilla",
        os.path.join(scripts_directory,"predict_from_epigenome.R"),
        f"--input_file={input_file}",
        f"--output_file={output_file}",
        f"--trained_model={model_path}"
    ]

    print(" ".join(com))

    subprocess.run(com)