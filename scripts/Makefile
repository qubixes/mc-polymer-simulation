.PHONY: clean
%.dvi: %.tex
	latex $< &> /dev/null
%.eps: %.dvi
	dvips -E -o $@ $< &> /dev/null
	
clean:
	rm -f *.{tex,eps,dvi,aux,log} *.tmp