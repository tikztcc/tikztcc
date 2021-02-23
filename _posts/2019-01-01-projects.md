---
title:  <h1 style="color:#000000;">Segunda imagem</h1>
#subtitle:
author: "João Gollner"
image: "img/imagensTCC/imgtcc2.png"
date:   2021-02-16 12:12:12
---


### Descrição

<p style="text-align: justify;">
Esta imagem foi retirada do livro Álgebra Linear do autor José Luiz Boldrini.
</p>

### Código para TIKZ no Latex

<p style="text-align: justify;">

\draw[->] (-2,0) -- (7,0) node[right]{$x_{1}$}; <br>
\draw[->] (0,-3) -- (0,6) node[right]{$x_{2}$}; <br>
\node[] at (-0.2,-0.2) {$0$}; <br>
\draw[domain=-1:4] plot (\x,{(-2*\x)+5}); <br>
\draw[domain=-3:7] plot (\x,{(-6+(\x))/3}); <br>
\draw (3,-1) node[right]{$(3,-1)$}; <br>
\draw (3,-1) node{$\bullet$};

</p>