import argparse
import os
from pyteomics import mzml
import json
from glob import glob

def summarize_scans(input_folder, output_file):
    output = []
    
    with open(output_file, 'w', encoding='utf-8') as json_file:
        for file in glob(input_folder + '/*.mzML'):            
            with mzml.read(file) as reader:
                for spectrum in reader:

                    output_dict = {}
                    output_dict['filename']         = os.path.basename(file)
                    output_dict['index']            = spectrum.get('index')
                    output_dict['id']               = spectrum.get('id')
                    output_dict['ms_level']         = spectrum.get('ms level')
                    if 'positive scan' in spectrum.keys():
                        output_dict['ion_mode']     = 'positive'
                    elif 'negative scan' in spectrum.keys():
                        output_dict['ion_mode']     = 'negative'
                    else:
                        output_dict['ion_mode']     = 'unknown'
                        
                    scan = None
                    try:
                        scan = spectrum["id"].replace("scanId=", "").split("scan=")[-1]
                    except Exception as _:
                        scan = spectrum["id"]
                    output_dict['scan'] = scan
                                       

                    mz = spectrum["m/z array"]
                    intensity = spectrum["intensity array"]
                    
                    # Discretize the m/z values to two decimal places, adding the intensities
                    mz_intensity = {}
                    for m, i in zip(mz, intensity):
                        mz_rounded = str(round(m, 2))
                        if mz_rounded in mz_intensity:
                            mz_intensity[mz_rounded] += i
                        else:
                            mz_intensity[mz_rounded] = i
                            
                    output_dict['m/z array'] = [str(x) for x in sorted(mz_intensity.keys())]
                    max_intensity = max(mz_intensity.values())
                    output_dict['intensity array'] = [str(round(mz_intensity[m]/max_intensity,4)) for m in output_dict['m/z array']]  # Normalize intensities
                    
                    output.append(output_dict)
                
        json.dump(output, json_file, indent=0)

def main():
    parser = argparse.ArgumentParser(description='Summarize m/z values for a folder of mzML files')
    parser.add_argument('--input_folder')
    parser.add_argument('--output_file')
    args = parser.parse_args()
    
    summarize_scans(args.input_folder, args.output_file)

if __name__ == "__main__":
    main()