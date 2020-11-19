// simple fragment shader
varying vec3 normal;
varying vec3 lightVec;
varying vec3 viewVec;

vec4 c1 = vec4(0.0,76,153,255);
vec4 c2 = vec4(0.0,128,255,255);
vec4 c3 = vec4(102,178,255,255);
vec4 c4 = vec4(204,229,255,255);
void main()
{
	vec3 norm = normalize(normal);
	vec3 L = normalize(lightVec);
	
	float NdotL = dot(norm,L);
	
	vec4 color;
	if(NdotL > 0.8)
	{
		color = c4/255;
	}
	else if(NdotL > 0.6)
	{
		color = c3/255;
	}
	else if(NdotL > 0.4)
	{
		color = c1/255;
	}
	else
	{
		color = c1/255;
	}
	
	gl_FragColor = color;
}
