import os,sys
import subprocess


def make_enpact_predictions(input_file, output_file, model_path):
    if os.path.exists(os.path.join(intermediates_dir,f"enpact_personalized_predictions_{ind}.txt")):
        print(ind," already exists")
        return

    com = [
        "Rscript",
        "--vanilla",
        os.path.join(script_directory,"predict_from_epigenome.R"),
        f"--input_file={input_file}",
        f"--output_file={output_file}",
        f"--trained_model={model_path}"
    ]

    print(" ".join(com))

    subprocess.run(com)