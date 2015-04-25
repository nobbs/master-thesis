.PHONY: main.pdf all clean

# The first rule in a Makefile is the one executed by default ("make"). It
# should always be the "all" rule, so that "make" and "make all" are identical.
all: main.pdf

view:
	open main.pdf

main.pdf: main.tex
	latexmk -silent -pdf -use-make main.tex

clean:
	latexmk -silent -CA

