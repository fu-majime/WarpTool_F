Texture2D tempbuffer : register(t0);
float4 psmain(float4 input : SV_Position) : SV_Target {
	float4 col = tempbuffer[int2(input.xy)];
	return col.a <= 1.0 ? col : col / col.a;
}
