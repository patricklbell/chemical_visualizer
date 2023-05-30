export const vec2 = (wsm, x, y) => {
  wsm.x = x;
  wsm.y = y;
  return wsm;
};

export const vec3 = (wsm, x, y, z) => {
  wsm.x = x;
  wsm.y = y;
  wsm.z = z;
  return wsm;
};

export const vec4 = (wsm, x, y, z, w) => {
  wsm.x = x;
  wsm.y = y;
  wsm.z = z;
  wsm.w = w;
  return wsm;
};
