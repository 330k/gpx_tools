function geoDistance(lat1, lon1, lat2, lon2){
  const a = 6378137.0;
  const f = 1 / 298.257223563;
  const epsilon = f * (2 - f) / ((1 - f) ** 2);
  
  const degree = Math.PI / 180.0;
  const sin = Math.sin;
  const cos = Math.cos;
  const sqrt = Math.sqrt;
  const tan = Math.tan;
  const atan = Math.atan;
  const abs = Math.abs;
  const asin = Math.asin;
  
  const radlat1 = lat1 * degree;
  const radlon1 = lon1 * degree;
  const radlat2 = lat2 * degree;
  const radlon2 = lon2 * degree;
  
  const l = radlon2 - radlon1;
  const lprime = ((l + 3 * (Math.PI)) % (2 * Math.PI)) - Math.PI;
  const L = Math.abs(l);
  const Lprime = Math.PI - L;
  
  const Delta = (l >= 0) ? (radlat2 - radlat1) : (radlat1 - radlat2);
  const Sigma = radlat1 + radlat2;
  
  const u1 = (l >= 0) ? atan((1 - f) * tan(radlat1)) : atan((1 - f) * tan(radlat2));
  const u2 = (l >= 0) ? atan((1 - f) * tan(radlat2)) : atan((1 - f) * tan(radlat1));
  
  const SigmaPrime = u1 + u2;
  const DeltaPrime = u2 - u1;
  const xi = cos(SigmaPrime / 2);
  const xiprime = sin(SigmaPrime / 2);
  const eta = sin(DeltaPrime / 2);
  const etaprime = cos(DeltaPrime / 2);
  
  const x = sin(u1) * sin(u2);
  const y = cos(u1) * cos(u2);
  
  const c = y * cos(L) + x;
  
  //ゾーンの計算
  let s;
  if(c >= 0){
    // Zone 1
    console.log("Zone 1");
    let theta = L * (1 + f * y);
    let g, h, sigma, J, K, gamma, Gamma, zeta, zetaprime, D, E, F, G;
    
    do {
      g = sqrt(eta * eta * (cos(theta/2)**2) + xi * xi * (sin(theta/2)**2));
      h = sqrt(etaprime * etaprime * (cos(theta/2)**2) + xiprime * xiprime * (sin(theta/2)**2));
      sigma = 2 * atan(g / h);
      J = 2 * g * h;
      K = h * h - g * g;
      gamma = y * sin(theta) / J;
      Gamma = 1 - gamma * gamma;
      zeta = Gamma * K - 2 * x;
      zetaprime = zeta + x;
      D = 0.25 * f * (1 - f) - (3 / 16) * f * f * Gamma;
      G = f * gamma * gamma * (1 - 2 * D * Gamma) + f * zetaprime * (sigma / J) * (1 - D * Gamma + 0.5 * f * gamma * gamma) + 0.25 * f * f * zeta * zetaprime;
      
      theta = theta - F / (1 - G);
    } while (Math.abs(F) > 1e-15);
    
    const n0 = epsilon * Gamma / ((sqrt(1 + epsilon * Gamma) + 1) ** 2);
    const A = (1 + n0) * (1 + 1.25 * n0 * n0);
    const B = epsilon * (1 - 3 * n0 * n0 / 8) / ((sqrt(1 + epsilon * Gamma) + 1) ** 2);
    
    s = (1 - f) * a * A * (sigma - B * J * (zeta - 0.25 * B * (K * (Gamma * Gamma - 2 * zeta * zeta) - (1 / 6) * B * zeta * (1 - 4 * K * K) * (3 * Gamma * Gamma - 4 * zeta * zeta))));
    
  }else{
    let theta;
    if(c >= -cos(3 * degree * cos(u1))){
      // Zone 2
      console.log("Zone 2");
      theta = Lprime;
    }else{
      // Zone 3
      const R = f * Math.PI * (cos(u1) ** 2) * (1 - 0.25 * f * (1 + f) * (sin(u1) ** 2) + (3 / 16) * f * f * (sin(u1) ** 4));
      const d1 = Lprime * cos(u1) - R;console.log(d1);
      const d2 = abs(SigmaPrime) + R;
      const q = Lprime / (f * Math.PI);
      const f1 = 0.25 * f * (1 + 0.5 * f);
      const gamma0 = q + f1 * q - f1 * (q ** 3);
      
      if(Sigma !== 0){
        // Zone 3a
        console.log("Zone 3a");
        const A0 = atan(d1 / d2);
        const B0 = asin(R / sqrt(d1 ** 2 + d2 ** 2));
        
        const psi = A0 + B0;
        const j = gamma0 / cos(u1);
        const k = (1 + f1) * abs(SigmaPrime) * (1 - f * y) / (f * Math.PI * y);
        const j1 = j / (1 + k / cos(psi));
        const psiprime1 = asin(j1);
        const psiprime2 = asin(cos(u1) / cos(u2) * j1);
        
        theta = 2 * atan(tan((psiprime1 + psiprime2) / 2) * sin(abs(SigmaPrime) / 2) / cos(DeltaPrime / 2));
        
      }else{
        if(d1 > 0){
          // Zone 3b1
          console.log("Zone 3b1");
          theta = Lprime;
          
        }else if(d1 === 0){
          // Zone 3b2
          console.log("Zone 3b2");
          const Gamma = sin(u1) ** 2;
          const n0 = epsilon * Gamma / ((sqrt(1 + epsilon * Gamma) + 1) ** 2);
          const A = (1 + n0) * (1 + 1.25 * n0 * n0);
          
          s = (1 - f) * a * A * Math.PI;
          
          return s;
        
        }else{
          // Zone 3b3
          console.log("Zone 3b3");
          
          let gamma = gamma0;
          let Gamma, D;
          
          for(let c = 0; c < 30; c++){
            Gamma = 1 - gamma * gamma;
            D = 0.25 * f * (1 + f) - (3 / 16) * f * f * Gamma;
            
            if(abs(q / (1 - D * Gamma) - gamma) > 1e-15){
              gamma = q - (1 - D * Gamma);
              console.log(gamma);
              
            }else{
              break;
            }
          }
          
          const n0 = epsilon * Gamma / ((sqrt(1 + epsilon * Gamma) + 1) ** 2);
          const A = (1 + n0) * (1 + 1.25 * n0 * n0);
          s = (1 - f) * a * A * Math.PI;
          return s;
        }
      }
    }
    
    let g, h, sigma, J, K, gamma, Gamma, zeta, zetaprime, D, E, F, G;
    
    do {
      g = sqrt(eta * eta * (sin(theta/2)**2) + xi * xi * (cos(theta/2)**2));
      h = sqrt(etaprime * etaprime * (sin(theta/2)**2) + xiprime * xiprime * (cos(theta/2)**2));
      sigma = 2 * atan(g / h);
      J = 2 * g * h;
      K = h * h - g * g;
      gamma = y * sin(theta) / J;
      Gamma = 1 - gamma * gamma;
      zeta = Gamma * K - 2 * x;
      zetaprime = zeta + x;
      D = 0.25 * f * (1 - f) - (3 / 16) * f * f * Gamma;
      G = f * gamma * gamma * (1 - 2 * D * Gamma) + f * zetaprime * (sigma / J) * (1 - D * Gamma + 0.5 * f * gamma * gamma) + 0.25 * f * f * zeta * zetaprime;
      
      theta = theta - F / (1 - G);
    } while (Math.abs(F) > 1e-15);
    
    const n0 = epsilon * Gamma / ((sqrt(1 + epsilon * Gamma) + 1) ** 2);
    const A = (1 + n0) * (1 + 1.25 * n0 * n0);
    const B = epsilon * (1 - 3 * n0 * n0 / 8) / ((sqrt(1 + epsilon * Gamma) + 1) ** 2);
    
    s = (1 - f) * a * A * (sigma - B * J * (zeta - 0.25 * B * (K * (Gamma * Gamma - 2 * zeta * zeta) - (1 / 6) * B * zeta * (1 - 4 * K * K) * (3 * Gamma * Gamma - 4 * zeta * zeta))));
    
    
  }
  
  return s;
}
