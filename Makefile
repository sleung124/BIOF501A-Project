# ================================================
# Makefile for removing files created from nextflow
# Adapted from Tony Liang's Makefile: https://github.com/tonyliang19/biof501-project/blob/main/Makefile
# ================================================
help:
	@echo "'make clean' clears out ALL nextflow results"
	@echo "  -- run when done work session"
	@echo "'make run_Rstudio' starts an rstudio instance"

clean:
	@echo "Removing intermediate files like logs and work/ ..."
	@rm -rf work/
	@rm -f .nextflow.log*
	@rm -rf .nextflow/
	@rm -rf results/
	@echo "done!"

build:
	@echo "Building docker image..."
	docker build -t sleung124/spatial-pipeline .

run_Rstudio:
	@echo "Running Rstudio..."
	docker run --rm -it -p 8787:8787 -e PASSWORD=123 -v `pwd -W`:/home/rstudio sleung124/visium_repo