.PHONY: paper clean sim sim-validate sim-verify all

# Default: build the paper
all: paper

# Build paper
paper:
	cd paper && latexmk -pdf main.tex

# Clean LaTeX auxiliary files
clean:
	cd paper && latexmk -c

# Run all simulations
sim: sim-validate sim-verify

# Monte Carlo validation of deterministic masking FIM recovery
sim-validate:
	Rscript research/validate_deterministic_masking.R

# Verify uniform masking is not FIM-optimal
sim-verify:
	Rscript research/verify_masking_optimality.R
