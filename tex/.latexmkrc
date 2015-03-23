$pdflatex = 'pdflatex -synctex=1 %O %S';

# add_cus_dep( 'glo', 'gls', 0, 'makeglossaries' );
# sub makeglossaries {
# system( "makeglossaries \"$_[0]\"" );
# }

 # Custom dependency and function for nomencl package 
 add_cus_dep( 'nlo', 'nls', 0, 'makenlo2nls' );
 sub makenlo2nls {
 system( "makeindex -s nomencl.ist -o \"$_[0].nls\" \"$_[0].nlo\"" );
 }
