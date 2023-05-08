run:
	nextflow run ./nf_workflow.nf -resume -c nextflow.config --input_spectra_folder=./data

run_hpcc:
	nextflow run ./nf_workflow.nf -resume -c nextflow_hpcc.config

run_docker:
	nextflow run ./nf_workflow.nf -resume -with-docker <CONTAINER NAME>