.PHONY: main.pdf all clean

# The first rule in a Makefile is the one executed by default ("make"). It
# should always be the "all" rule, so that "make" and "make all" are identical.
all: main.pdf

view:
	latexmk	-silent late-pv -pdf main.pdf

main.pdf: main.tex
	latexmk -silent -pdf main.tex

continuous:
	latexmk -pvc -pdf -silent -f main.tex

clean:
	latexmk -silent -CA

clean-biber:
	rm -rf `biber --cache`

clean-all: clean clean-biber

