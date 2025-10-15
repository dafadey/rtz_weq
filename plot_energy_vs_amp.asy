import graph;

void read_point(real[][] data, string filename) {
  file inf = input(filename).line();
  string[] raw = inf;
  //write(raw);
  string[] sinput = split(raw[4], '\t');
  real[] input;
  for(string v : sinput)
    input.push((real) v);
  data.push(input);
}

pair data[];

string number(int i)
{
  string num;
  string _num = string(i);
  for(int j=0; j<7-length(_num);++j)
    num += '0';
  num += _num;
  return num;
}

void add_point(real amp, string name){
  real[][] d;
  read_point(d, name+"/"+"field_energy"+number(95)+".dat");
  //data.push((amp, d[0][2]/d[0][0]));
  data.push((d[0][0], d[0][2]));
}

add_point(0.00011714, "w_harm_ampW0=0.00011744_ratio2w=0_/pics/res");
add_point(0.0002, "w_harm_ampW0=0.0002_ratio2w=0_/pics/res");
add_point(0.0003, "w_harm_ampW0=0.0003_ratio2w=0_/pics/res");
add_point(0.0004, "w_harm_ampW0=0.0004_ratio2w=0_/pics/res");
add_point(0.0005, "w_harm_ampW0=0.0005_ratio2w=0_/pics/res");
add_point(0.0006, "w_harm_ampW0=0.0006_ratio2w=0_/pics/res");
add_point(0.0007, "w_harm_ampW0=0.0007_ratio2w=0_/pics/res");
add_point(0.0008, "w_harm_ampW0=0.0008_ratio2w=0_/pics/res");
add_point(0.0009,"w_harm_ampW0=0.0009_ratio2w=0_/pics/res");
add_point(0.0010, "w_harm_ampW0=0.0010_ratio2w=0_/pics/res");
add_point(0.0011, "w_harm_ampW0=0.0011_ratio2w=0_/pics/res");
add_point(0.0011744, "w_harm_ampW0=0.0011744_ratio2w=0_/pics/res");

guide g;
guide gd;
for(pair item : data) {
  write(item);
  g=g--(log10(item.x), log10(item.y));
}

for(int i=1;i<data.length;++i)
  
  gd=gd--((log10(data[i].x)+log10(data[i-1].x))*.5, (log10(data[i].y)-log10(data[i-1].y))/(log10(data[i].x)-log10(data[i-1].x)));

picture pic;
scale(pic, Log, Log);
draw(pic, g);

xaxis(pic, "$W_{\omega}$", BottomTop, LeftTicks);
yaxis(pic, "$W_{3\omega}$", LeftRight, RightTicks);
size(pic, 8cm, 5cm, point(pic, SW), point(pic, NE));
pic = shift(-point(pic,SW))*pic;
add(pic.fit());

picture picd;
scale(picd, Log, Linear);
draw(picd, gd);

xaxis(picd, "$W_{\omega}$", BottomTop, LeftTicks);
yaxis(picd, "slope", LeftRight, RightTicks);
size(picd, 8cm, 5cm, point(picd, SW), point(picd, NE));
picd = shift(-point(picd,SW))*picd;
add(picd.fit(),(0,-7cm));

