#include <stdio.h>
#include <assert.h>
#define TRUE 1
#define FALSE 0
void top(uint32_t x[0], uint32_t x[1], uint32_t x[2], uint32_t x[3], uint32_t y[0], uint32_t y[1], uint32_t y[2], uint32_t y[3], uint32_t *z[0], uint32_t *z[1], uint32_t *z[2], uint32_t *z[3])
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
  n17 = x[3] & y[3];
  n16 = x[2] & y[2];
  n21 = n17 ^ n16;
  n10 = x[1] & y[1];
  n9 = x[0] & y[0];
  n14 = n10 ^ n9;
  n35 = n21 ^ n14;
  n18 = x[3] ^ x[2];
  n19 = y[3] ^ y[2];
  n20 = n18 & n19;
  n22 = n20 ^ n16;
  n11 = x[1] ^ x[0];
  n12 = y[1] ^ y[0];
  n13 = n11 & n12;
  n15 = n13 ^ n9;
  n36 = n22 ^ n15;
  n25 = x[3] ^ x[1];
  n26 = y[3] ^ y[1];
  n29 = n25 & n26;
  n23 = x[2] ^ x[0];
  n24 = y[2] ^ y[0];
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
  *z[0] = n35;
  *z[1] = n36;
  *z[2] = n38;
  *z[3] = n40;
}