 
/************* HSI to RGB conversion routine ****************/
// R,G,B values are from 0 to 1
// H = [0,360], S = [0,1], V = [0,1]

void HSVtoRGB( float RGB[3], float h, float s, float v )
{

  int i;
  double f;


  h = h/120.0;
  i = floor( h);
  f = h - (double) i;

  if( s == 0 ) {
          // achromatic (grey)
    RGB[0] = RGB[1] = RGB[2] = v;
    return;
  }
//printf("%d ", i);

  switch (i)  {
    case 0: 
      RGB[0] = v * (1.0-s + s * (1.0-f));
      RGB[1] = v * (1.0-s + s * f);
      RGB[2] = v * (1.0-s);
      break;
    case 1: 
      RGB[1] = v * (1.0-s + s * (1.0-f));
      RGB[2] = v * (1.0-s + s * f);
      RGB[0] = v * (1.0-s);
      break;
    default :
      RGB[2] = v * (1.0-s + s * (1.0-f));
      RGB[0] = v * (1.0-s + s * f);
      RGB[1] = v * (1.0-s);
      break;
  }

}

