module top( x0 , x1 , x2 , x3 , x4 , x5 , x6 , x7 , y0 , y1 , y2 , y3 );
  input x0 , x1 , x2 , x3 , x4 , x5 , x6 , x7 ;
  output y0 , y1 , y2 , y3 ;
  wire n9 , n10 , n11 , n12 , n13 , n14 , n15 , n16 , n17 , n18 , n19 , n20 , n21 , n22 , n23 , n24 , n25 , n26 , n27 , n28 , n29 , n30 , n31 , n32 , n33 , n34 , n35 , n36 , n37 , n38 , n39 , n40 ;
  assign n17 = x3 & x7 ;
  assign n16 = x2 & x6 ;
  assign n21 = n17 ^ n16 ;
  assign n10 = x1 & x5 ;
  assign n9 = x0 & x4 ;
  assign n14 = n10 ^ n9 ;
  assign n35 = n21 ^ n14 ;
  assign n18 = x3 ^ x2 ;
  assign n19 = x7 ^ x6 ;
  assign n20 = n18 & n19 ;
  assign n22 = n20 ^ n16 ;
  assign n11 = x1 ^ x0 ;
  assign n12 = x5 ^ x4 ;
  assign n13 = n11 & n12 ;
  assign n15 = n13 ^ n9 ;
  assign n36 = n22 ^ n15 ;
  assign n25 = x3 ^ x1 ;
  assign n26 = x7 ^ x5 ;
  assign n29 = n25 & n26 ;
  assign n23 = x2 ^ x0 ;
  assign n24 = x6 ^ x4 ;
  assign n28 = n23 & n24 ;
  assign n33 = n29 ^ n28 ;
  assign n37 = n35 ^ n33 ;
  assign n38 = n37 ^ n22 ;
  assign n30 = n25 ^ n23 ;
  assign n31 = n26 ^ n24 ;
  assign n32 = n30 & n31 ;
  assign n34 = n32 ^ n28 ;
  assign n39 = n36 ^ n34 ;
  assign n27 = n22 ^ n21 ;
  assign n40 = n39 ^ n27 ;
  assign y0 = n35 ;
  assign y1 = n36 ;
  assign y2 = n38 ;
  assign y3 = n40 ;
endmodule
