$pdflatex = 'pdflatex -synctex=1 %O %S';
$clean_ext = 'synctex.gz synctex.gz(busy) run.xml tex.bak bbl bcf fdb_latexmk run tdo %R-blx.bib bcf nlo thm sta paux lox pdfsync tdo syi syg slg';
$pdf_previewer = 'open -a /Applications/Skim.app';
# $aux_dir = 'build';

# Custom dependency for glossary
# add_cus_dep( 'glo', 'gls', 0, 'makeglossaries' );
# sub makeglossaries {
# system( "makeglossaries \"$_[0]\"" );
# }

# add_cus_dep( 'syg', 'syi', 0, 'makeglossaries' );
# sub makeglossaries {
# system( "makeglossaries \"$_[0]\"" );
# }


# Custom dependency and function for nomencl package
# add_cus_dep( 'nlo', 'nls', 0, 'makenlo2nls' );
# sub makenlo2nls {
# system( "makeindex -s nomencl.ist -o \"$_[0].nls\" \"$_[0].nlo\"" );
# }
