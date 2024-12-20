#include <stdio.h>
#include <assert.h>
#define TRUE 1
#define FALSE 0
void top(uint32_t x0, uint32_t x1, uint32_t x2, uint32_t x3, uint32_t x4, uint32_t x5, uint32_t x6, uint32_t x7, uint32_t *y0, uint32_t *y1, uint32_t *y2, uint32_t *y3)
{
  uint32_t n9;
  uint32_t n10;
  uint32_t n11;
  uint32_t n12;
  uint32_t n13;
  uint32_t n14;
  uint32_t n15;
  uint32_t n16;
  uint32_t n17;
  uint32_t n18;
  uint32_t n19;
  uint32_t n20;
  uint32_t n21;
  uint32_t n22;
  uint32_t n23;
  uint32_t n24;
  uint32_t n25;
  uint32_t n26;
  uint32_t n27;
  uint32_t n28;
  uint32_t n29;
  uint32_t n30;
  uint32_t n31;
  uint32_t n32;
  uint32_t n33;
  uint32_t n34;
  uint32_t n35;
  uint32_t n36;
  uint32_t n37;
  uint32_t n38;
  uint32_t n39;
  uint32_t n40;
  n17 = x3 & x7;
  n16 = x2 & x6;
  n21 = n17 ^ n16;
  n10 = x1 & x5;
  n9 = x0 & x4;
  n14 = n10 ^ n9;
  n35 = n21 ^ n14;
  n18 = x3 ^ x2;
  n19 = x7 ^ x6;
  n20 = n18 & n19;
  n22 = n20 ^ n16;
  n11 = x1 ^ x0;
  n12 = x5 ^ x4;
  n13 = n11 & n12;
  n15 = n13 ^ n9;
  n36 = n22 ^ n15;
  n25 = x3 ^ x1;
  n26 = x7 ^ x5;
  n29 = n25 & n26;
  n23 = x2 ^ x0;
  n24 = x6 ^ x4;
  n28 = n23 & n24;
  n33 = n29 ^ n28;
  n37 = n35 ^ n33;
  n38 = n37 ^ n22;
  n30 = n25 ^ n23;
  n31 = n26 ^ n24;
  n32 = n30 & n31;
  n34 = n32 ^ n28;
  n39 = n36 ^ n34;
  n27 = n22 ^ n21;
  n40 = n39 ^ n27;
  *y0 = n35;
  *y1 = n36;
  *y2 = n38;
  *y3 = n40;
}
