CPPFLAGS = -I../dmrg/Engine -I../dmrg/Models/Heisenberg -I../dmrg/Models/HubbardOneBand -I../dmrg/Models/FeAsModel \
	   -I../../PsimagLite/src -I../../PsimagLite/src/Geometry
CXX = g++

all: manual.pdf

FILES = $(shell find ../dmrg -iname "*.h" -or -iname "*.cpp")
manual.tex: manual.ptex ../README.md  $(FILES)
	find ../dmrg -iname "*.h" -or -iname "*.cpp" |\
              perl ../PsimagLite/scripts/doc.pl manual.ptex

manual.pdf: manual.tex ../README.md ../dmrg/*h ../dmrg/Engine/*h ../dmrg/*cpp
	pdflatex manual.tex
	bibtex manual.aux
	pdflatex manual.tex
	pdflatex manual.tex

clean:
	rm -f ../README.tex
	rm -f manual.tex manual.pdf


