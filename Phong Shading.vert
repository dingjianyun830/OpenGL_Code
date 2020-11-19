// simple vertex shader

varying vec3 normal;
varying vec4 vert;

void main()
{
	gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
	vert = gl_ModelViewMatrix * gl_Vertex;
	normal = gl_NormalMatrix *gl_Normal;
}
