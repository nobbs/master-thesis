$pdflatex = 'pdflatex -synctex=1 %O %S';

add_cus_dep( 'glo', 'gls', 0, 'makeglossaries' );
sub makeglossaries {
   system( "makeglossaries \"$_[0]\"" );
}
