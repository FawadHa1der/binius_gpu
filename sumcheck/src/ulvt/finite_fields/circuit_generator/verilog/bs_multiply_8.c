#include <stdio.h>
#include <assert.h>
#define TRUE 1
#define FALSE 0
void top(_Bool x0, _Bool x1, _Bool x2, _Bool x3, _Bool x4, _Bool x5, _Bool x6, _Bool x7, _Bool x8, _Bool x9, _Bool x10, _Bool x11, _Bool x12, _Bool x13, _Bool x14, _Bool x15, _Bool *y0, _Bool *y1, _Bool *y2, _Bool *y3, _Bool *y4, _Bool *y5, _Bool *y6, _Bool *y7)
{
  _Bool n17;
  _Bool n18;
  _Bool n19;
  _Bool n20;
  _Bool n21;
  _Bool n22;
  _Bool n23;
  _Bool n24;
  _Bool n25;
  _Bool n26;
  _Bool n27;
  _Bool n28;
  _Bool n29;
  _Bool n30;
  _Bool n31;
  _Bool n32;
  _Bool n33;
  _Bool n34;
  _Bool n35;
  _Bool n36;
  _Bool n37;
  _Bool n38;
  _Bool n39;
  _Bool n40;
  _Bool n41;
  _Bool n42;
  _Bool n43;
  _Bool n44;
  _Bool n45;
  _Bool n46;
  _Bool n47;
  _Bool n48;
  _Bool n49;
  _Bool n50;
  _Bool n51;
  _Bool n52;
  _Bool n53;
  _Bool n54;
  _Bool n55;
  _Bool n56;
  _Bool n57;
  _Bool n58;
  _Bool n59;
  _Bool n60;
  _Bool n61;
  _Bool n62;
  _Bool n63;
  _Bool n64;
  _Bool n65;
  _Bool n66;
  _Bool n67;
  _Bool n68;
  _Bool n69;
  _Bool n70;
  _Bool n71;
  _Bool n72;
  _Bool n73;
  _Bool n74;
  _Bool n75;
  _Bool n76;
  _Bool n77;
  _Bool n78;
  _Bool n79;
  _Bool n80;
  _Bool n81;
  _Bool n82;
  _Bool n83;
  _Bool n84;
  _Bool n85;
  _Bool n86;
  _Bool n87;
  _Bool n88;
  _Bool n89;
  _Bool n90;
  _Bool n91;
  _Bool n92;
  _Bool n93;
  _Bool n94;
  _Bool n95;
  _Bool n96;
  _Bool n97;
  _Bool n98;
  _Bool n99;
  _Bool n100;
  _Bool n101;
  _Bool n102;
  _Bool n103;
  _Bool n104;
  _Bool n105;
  _Bool n106;
  _Bool n107;
  _Bool n108;
  _Bool n109;
  _Bool n110;
  _Bool n111;
  _Bool n112;
  _Bool n113;
  _Bool n114;
  _Bool n115;
  _Bool n116;
  _Bool n117;
  _Bool n118;
  _Bool n119;
  _Bool n120;
  _Bool n121;
  _Bool n122;
  _Bool n123;
  _Bool n124;
  _Bool n125;
  _Bool n126;
  _Bool n127;
  _Bool n128;
  _Bool n129;
  _Bool n130;
  _Bool n131;
  _Bool n132;
  _Bool n133;
  _Bool n134;
  _Bool n135;
  n57 = x7 && x15;
  n56 = x6 && x14;
  n61 = n57 ^ n56;
  n50 = x5 && x13;
  n49 = x4 && x12;
  n54 = n50 ^ n49;
  n75 = n61 ^ n54;
  n25 = x3 && x11;
  n24 = x2 && x10;
  n29 = n25 ^ n24;
  n18 = x1 && x9;
  n17 = x0 && x8;
  n22 = n18 ^ n17;
  n43 = n29 ^ n22;
  n124 = n75 ^ n43;
  n58 = x7 ^ x6;
  n59 = x15 ^ x14;
  n60 = n58 && n59;
  n62 = n60 ^ n56;
  n51 = x5 ^ x4;
  n52 = x13 ^ x12;
  n53 = n51 && n52;
  n55 = n53 ^ n49;
  n76 = n62 ^ n55;
  n26 = x3 ^ x2;
  n27 = x11 ^ x10;
  n28 = n26 && n27;
  n30 = n28 ^ n24;
  n19 = x1 ^ x0;
  n20 = x9 ^ x8;
  n21 = n19 && n20;
  n23 = n21 ^ n17;
  n44 = n30 ^ n23;
  n125 = n76 ^ n44;
  n65 = x7 ^ x5;
  n66 = x15 ^ x13;
  n69 = n65 && n66;
  n63 = x6 ^ x4;
  n64 = x14 ^ x12;
  n68 = n63 && n64;
  n73 = n69 ^ n68;
  n77 = n75 ^ n73;
  n78 = n77 ^ n62;
  n33 = x3 ^ x1;
  n34 = x11 ^ x9;
  n37 = n33 && n34;
  n31 = x2 ^ x0;
  n32 = x10 ^ x8;
  n36 = n31 && n32;
  n41 = n37 ^ n36;
  n45 = n43 ^ n41;
  n46 = n45 ^ n30;
  n126 = n78 ^ n46;
  n70 = n65 ^ n63;
  n71 = n66 ^ n64;
  n72 = n70 && n71;
  n74 = n72 ^ n68;
  n79 = n76 ^ n74;
  n67 = n62 ^ n61;
  n80 = n79 ^ n67;
  n38 = n33 ^ n31;
  n39 = n34 ^ n32;
  n40 = n38 && n39;
  n42 = n40 ^ n36;
  n47 = n44 ^ n42;
  n35 = n30 ^ n29;
  n48 = n47 ^ n35;
  n127 = n80 ^ n48;
  n87 = x7 ^ x3;
  n88 = x15 ^ x11;
  n100 = n87 && n88;
  n85 = x6 ^ x2;
  n86 = x14 ^ x10;
  n99 = n85 && n86;
  n104 = n100 ^ n99;
  n83 = x5 ^ x1;
  n84 = x13 ^ x9;
  n93 = n83 && n84;
  n81 = x4 ^ x0;
  n82 = x12 ^ x8;
  n92 = n81 && n82;
  n97 = n93 ^ n92;
  n118 = n104 ^ n97;
  n128 = n124 ^ n118;
  n129 = n128 ^ n78;
  n101 = n87 ^ n85;
  n102 = n88 ^ n86;
  n103 = n101 && n102;
  n105 = n103 ^ n99;
  n94 = n83 ^ n81;
  n95 = n84 ^ n82;
  n96 = n94 && n95;
  n98 = n96 ^ n92;
  n119 = n105 ^ n98;
  n130 = n125 ^ n119;
  n131 = n130 ^ n80;
  n108 = n87 ^ n83;
  n109 = n88 ^ n84;
  n112 = n108 && n109;
  n106 = n85 ^ n81;
  n107 = n86 ^ n82;
  n111 = n106 && n107;
  n116 = n112 ^ n111;
  n120 = n118 ^ n116;
  n121 = n120 ^ n105;
  n132 = n126 ^ n121;
  n90 = n80 ^ n75;
  n133 = n132 ^ n90;
  n113 = n108 ^ n106;
  n114 = n109 ^ n107;
  n115 = n113 && n114;
  n117 = n115 ^ n111;
  n122 = n119 ^ n117;
  n110 = n105 ^ n104;
  n123 = n122 ^ n110;
  n134 = n127 ^ n123;
  n89 = n80 ^ n78;
  n91 = n89 ^ n76;
  n135 = n134 ^ n91;
  *y0 = n124;
  *y1 = n125;
  *y2 = n126;
  *y3 = n127;
  *y4 = n129;
  *y5 = n131;
  *y6 = n133;
  *y7 = n135;
}
void main() {
_Bool x0;
_Bool x1;
_Bool x2;
_Bool x3;
_Bool x4;
_Bool x5;
_Bool x6;
_Bool x7;
_Bool x8;
_Bool x9;
_Bool x10;
_Bool x11;
_Bool x12;
_Bool x13;
_Bool x14;
_Bool x15;
_Bool y0;
_Bool y1;
_Bool y2;
_Bool y3;
_Bool y4;
_Bool y5;
_Bool y6;
_Bool y7;
top(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, &y0, &y1, &y2, &y3, &y4, &y5, &y6, &y7);
}
