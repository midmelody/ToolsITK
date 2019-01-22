//Compute ellipsoid radius along neighborhood directions
float elliRadius(VECTOR vec1, VECTOR vec2, VECTOR vec3, int x, int y, int z)
{   
  VECTOR t,v; //temp variable
  float a, b, c; //length of three axis
  float sinTheta, cosTheta, sinPhi, cosPhi; //sin and cos
  float sin2Theta, cos2Theta, sin2Phi, cos2Phi; //square of sin and cos
  float elliLength;
  float norm;
  a = vec1.length;
  b = vec2.length;
  c = vec3.length;

  //adjacency vec
  t.x = -float(x);
  t.y = -float(y);
  t.z = -float(z);
  norm = sqrt(t.x*t.x+t.y*t.y+t.z*t.z);
  t.x = t.x/norm;
  t.y = t.y/norm;
  t.z = t.z/norm;

  sinTheta = t.x*vec3.x+t.y*vec3.y+t.z*vec3.z;
  sin2Theta = sinTheta*sinTheta;
  cos2Theta = 1-sin2Theta;

  v.x = t.x-vec3.x*sinTheta;
  v.y = t.y-vec3.y*sinTheta;
  v.z = t.z-vec3.z*sinTheta;
	
  if(sqrt(v.x*v.x+v.y*v.y+v.z*v.z)!=0)
    {
      cosPhi = (v.x*vec1.x+v.y*vec1.y+v.z*vec1.z)/sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
      cos2Phi = cosPhi*cosPhi;
      sin2Phi = 1-cos2Phi;
    }
  else
    {
      sin2Theta = 1;
      cos2Theta = 0;
      cos2Phi = 1;
      sin2Phi = 0;
    }
  elliLength = (a*a*b*b*c*c)/(a*a*b*b*sin2Theta+c*c*cos2Theta*(a*a*sin2Phi+b*b*cos2Phi));
  elliLength = sqrt(elliLength); 
  return(elliLength);
}
