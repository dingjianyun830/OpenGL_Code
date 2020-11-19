// simple vertex shader

attribute vec4 attr_c;// per vertex color

varying vec4 color;// pass into fragment

void main()
{
	gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
	vec4 vert = gl_ModelViewMatrix * gl_Vertex;
	vec3 normal = gl_NormalMatrix *gl_Normal;
	vec3 lightVec = vec3(gl_LightSource[0].position - vert);
	vec3 viewVec = vec3(-vert);

	vec3 norm = normalize(normal);
	vec3 L = normalize(lightVec);
	vec3 view = normalize(viewVec);

	vec3 refl = normalize(reflect(L,norm));
	float NdotL = dot(norm,L);
	float VdotR = dot(view,refl);
	
	vec4 ambient  = gl_LightSource[0].ambient * gl_FrontMaterial.ambient;
	vec4 diffuse  = gl_LightSource[0].diffuse * gl_FrontMaterial.diffuse * NdotL + 0.4;
    vec4 specular = gl_LightSource[0].specular * gl_FrontMaterial.specular * pow(max(0.0,VdotR),gl_FrontMaterial.shininess);
	vec4 Itotal = ambient + diffuse + specular;	

	color = attr_c*Itotal;
}
