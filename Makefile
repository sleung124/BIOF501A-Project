# ================================================
# Makefile for removing files created from nextflow
# Adapted from Tony Liang's Makefile: https://github.com/tonyliang19/biof501-project/blob/main/Makefile
# ================================================
help:
	@echo "'make clean' clears out ALL nextflow results"
	@echo "  -- run when done work session"

clean:
	@echo "Removing intermediate files like logs and work/ ..."
	@rm -rf work/
	@rm -f .nextflow.log*
	@rm -rf .nextflow/
	@rm -rf results/
	@echo "done!"
