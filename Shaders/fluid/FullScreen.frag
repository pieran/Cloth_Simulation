#version 150 core

uniform sampler2D texNormals;
uniform sampler2D texAbsorbtion;
uniform sampler2D texDepth;
uniform samplerCube texReflect;

uniform vec3  ambientColour;
uniform vec3  invLightDir;			//Directional Light
uniform vec3  cameraPos;
uniform float specularIntensity = 64.0f;
uniform float baseReflect = 0.1f;
//uniform float baseRefract = 0.8f;

uniform vec3 absorbtionExp = vec3(5.5f, 3.2f, 0.1f);
uniform mat4 invProjViewMatrix;


in Vertex	{
	vec2 texCoords;
} IN;

out vec4 gl_FragColor;


void main(void)	{
	float ndc_depth = texture(texDepth, IN.texCoords).x;
	if (ndc_depth > 0.99f)
	{
		discard;
		return;
	}
	
	//Set Fragment Depth
	gl_FragDepth = ndc_depth;
	
	
	//Read Texture Specific Data
	vec3 normal = texture(texNormals, IN.texCoords).rgb * 2.0f - 1.0f;
	float absorbtion = texture(texAbsorbtion, IN.texCoords).x;
	
	
	//Compute World-Space Position
	vec4 clip_pos = vec4(IN.texCoords.x, IN.texCoords.y, ndc_depth, 1.0f) * 2.0f - 1.0f;
	vec4 ws_pos = invProjViewMatrix * clip_pos;
	ws_pos.xyz /= ws_pos.w;
	
	
	
	//Lighting Calculations (Diffuse + Specular Factors)
	vec3 viewDir 		= normalize(cameraPos - ws_pos.xyz );
	vec3 halfDir 		= normalize(invLightDir + viewDir );
	float rFactor       = max(0.0, dot(halfDir , normal ));
	
	float dFactor       = max(0.0, dot(invLightDir , normal )) ;
    float sFactor       = pow(rFactor , specularIntensity );
	   
	//Lighting Calculations (Fresnel) 
	float theta = clamp(dot(viewDir, normal), 0.0f, 1.0f);
	float fFactor = baseReflect + (1.0f - baseReflect) * pow(1.0f - theta, 2.0f);  
	   

	//Calculate Base Color (from absorbtion)
	vec3 diffColour = vec3(1.0f);
	diffColour.r = exp(-absorbtionExp.r * absorbtion);
	diffColour.g = exp(-absorbtionExp.g * absorbtion);
	diffColour.b = exp(-absorbtionExp.b * absorbtion);
	
	
	float ratio = 1.0 /1.3333; //refractive index
	vec3 refractCoords = refract(-viewDir, normal, ratio);
	vec3 viewDir2 		= -normalize(ws_pos.xyz - cameraPos);
	vec3 reflColour = texture(texReflect, normalize(reflect(-viewDir, normal))).rgb;
	vec3 refrColour = texture(texReflect, refractCoords).rgb;// - absorbtion * normal * 0.1f).rgb;
	
	float baseRefract = 1.0f - clamp(exp(-2.15f * (1.0f - absorbtion * 0.1f)), 0.0f, 1.0f);	
	diffColour *= (1.0f - baseRefract);
	diffColour += refrColour * baseRefract;
	diffColour *= (1.0f - fFactor);
	diffColour += reflColour * fFactor;
//diffColour = refrColour;
	
	//Colour Computations
	vec3 specColour = min(diffColour + vec3(0.5f), vec3(1)); //Quick hack to approximate specular colour of an object, assuming the light colour is white
    vec3 col = ambientColour * diffColour 
			 + diffColour * dFactor
			 + specColour * sFactor;
		
		
	//Output Final Fragment Colour
	gl_FragColor.rgb 	= col;
	gl_FragColor.a 		= 1.0f - clamp(exp(-5.5f * absorbtion), 0.0f, 1.0f);
}