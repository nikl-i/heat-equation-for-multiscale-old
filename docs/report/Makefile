main := src/main.tex

world: pdf
	
pdf:
	pdflatex $(main)
dvi:
	latex $(main)
clean:
	rm -r auto/ *.pdf *.log *.dvi *.aux *.out *~ "#main.tex#"
