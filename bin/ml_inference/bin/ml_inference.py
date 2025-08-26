import argparse
from pathlib import Path
import pandas as pd 
from torchvision import transforms
import sys
import onnxruntime as ort
import numpy as np
from tqdm import tqdm
from torch import Tensor
from custom_transforms import SquareRootTransform, SelectTopKPeaks, BinarizeIntensity, NormalizeIntensity, PadToLength

trans = transforms.Compose([
    SquareRootTransform(),
    SelectTopKPeaks(150),
    BinarizeIntensity(), # *************
    NormalizeIntensity(),
    PadToLength(150, padding_value=-1.0),
])

def run_inference(ml_data_directory: Path, output_file: Path, session: ort.InferenceSession):

    outputs = {}

    for file in tqdm(list(ml_data_directory.glob('*.npy'))):
        if not file.is_file():
            continue
        
        # Load the data
        data = np.load(file)
        db_id = file.stem

        transformed = trans(Tensor(data).T).unsqueeze(0).numpy().astype(np.float32)

        outputs[db_id] = session.run(None, {session.get_inputs()[0].name: transformed})[0].squeeze().tolist()


    print(outputs)
    # Convert outputs to DataFrame
    df = pd.DataFrame.from_dict(outputs, orient='index')
    df.reset_index(inplace=True)
    df.rename(columns={'index': 'database_id'}, inplace=True)
    df.to_feather(output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ML inference on preprocessed data")
    parser.add_argument('--ml_data_directory', type=str, required=True, help='Directory containing preprocessed data')
    parser.add_argument('--output_file', type=str, required=True, help='Output file to save the inference results')
    parser.add_argument('--onnx_model', type=str, required=True, help='Path to the ONNX model file')
    args = parser.parse_args()

    ml_data_directory = Path(args.ml_data_directory)
    output_file = Path(args.output_file)
    onnx_model = Path(args.onnx_model)

    if not ml_data_directory.exists():
        raise ValueError(f"ML data directory does not exist: {ml_data_directory}")
    if not onnx_model.exists():
        raise ValueError(f"ONNX model file does not exist: {onnx_model}")
    if not output_file.parent.exists():
        output_file.parent.mkdir(parents=True, exist_ok=True)

    session = ort.InferenceSession(onnx_model, providers=["CPUExecutionProvider"])

    # Get input and output names
    input_name = session.get_inputs()[0].name
    output_name = session.get_outputs()[0].name

    run_inference(ml_data_directory, output_file, session)