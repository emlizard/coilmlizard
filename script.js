// Elliptic integrals via AGM algorithm
    function ellipke(m) {
      let a = 1.0;
      let b = Math.sqrt(1 - m);
      let sumE = 1 - m / 2;
      let iter = 0;
      while (Math.abs(a - b) > 1e-12 && iter < 50) {
        let an = (a + b) / 2.0;
        let bn = Math.sqrt(a * b);
        let cn = (a - b) / 2.0;
        sumE -= Math.pow(2, iter) * cn * cn;
        a = an; b = bn;
        iter++;
      }
      let K = Math.PI / (2 * a);
      let E = K * sumE;
      return [K, E];
    }

    // f24 integrand evaluation over array of p
    function f24Array(pArr, alpha, beta, gamma, delta, a, b, c) {
      const results = new Array(pArr.length);
      const h = alpha, e = beta, g = gamma, d = delta;
      const a2 = a * a, b2 = b * b, c2 = c * c;
      const l2 = a2 + c2, l = Math.sqrt(l2);
      const L2 = l2 + b2, L = Math.sqrt(L2);
      const l2L2 = l2 * L2, lL = l * L;

      for (let idx = 0; idx < pArr.length; idx++) {
        const p = pArr[idx];
        const sp = Math.sin(p), cp = Math.cos(p);
        let V, p1, p2, p3, p4, p5;

        if (l === 0) {
          p1 = 0;
          p2 = -g * Math.sign(b);
          p3 = 0;
          p4 = -e * Math.sign(b);
          p5 = d;
          V = Math.sqrt(e*e + g*g + h*h*cp*cp - 2*h*e*Math.sign(b)*cp);
        } else {
          p1 = (g * c) / l;
          p2 = - (e * l2 + g * a * b) / (l * L);
          p3 = (h * c) / L;
          p4 = (g * l2 - e * a * b - d * b * c) / L;
          p5 = (d * a - e * c) / l;
          const term1 = e*e + g*g;
          const term2 = h*h * (
            (1 - (b2 * c2) / (l2L2)) * cp*cp
            + (c2 / l2) * sp*sp
            + (a * b * c) / (l2 * L) * Math.sin(2*p)
          );
          const term3 = - (2 * h / (lL)) * (e * a * b - g * l2) * cp;
          const term4 = - (2 * h * e * c / l) * sp;
          V = Math.sqrt(term1 + term2 + term3 + term4);
        }

        const A = (1 + e*e + g*g + h*h + d*d) + 2*h*(p4 * cp + p5 * sp);
        const m = 4 * V / (A + 2 * V);
        const k = Math.sqrt(m);
        const [Kval, Eval] = ellipke(m);
        const PSI = (2 - m) * Kval - 2 * Eval;
        results[idx] = (p1 * cp + p2 * sp + p3) * (PSI / (k * Math.pow(V, 1.5)));
      }
      return results;
    }

    // Babic_24 mutual inductance calculation
    function Babic_24(Rp, Rs, pc, n, tol = 1e-13) {
      const xc = pc[0] / Rp, yc = pc[1] / Rp, zc = pc[2] / Rp;
      const a = n[0], b = n[1], c = n[2];
      const alpha = Rs / Rp, beta = xc, gamma = yc, delta = zc;

      const decdigs = Math.abs(Math.floor(Math.log10(tol)));
      const nPts = Math.pow(2, decdigs - 1) + 1;
      const pArr = new Array(nPts);
      for (let i = 0; i < nPts; i++) {
        pArr[i] = (2 * Math.PI) * i / (nPts - 1);
      }
      const romall = f24Array(pArr, alpha, beta, gamma, delta, a, b, c);

      let romPrev = new Array(decdigs).fill(0);
      let romCurr = new Array(decdigs).fill(0);

      let h = 2 * Math.PI;
      romPrev[0] = h * (romall[0] + romall[nPts - 1]) / 2;

      for (let i = 2; i <= decdigs; i++) {
        const step = Math.pow(2, decdigs - i + 1);
        let sumNew = 0;
        for (let j = step / 2; j < nPts - 1; j += step) {
          sumNew += romall[j];
        }
        romCurr[0] = (romPrev[0] + h * sumNew) / 2;
        for (let k = 1; k < i; k++) {
          const factor = Math.pow(4, k);
          romCurr[k] = (factor * romCurr[k - 1] - romPrev[k - 1]) / (factor - 1);
        }
        for (let k = 0; k < i; k++) {
          romPrev[k] = romCurr[k];
        }
        h /= 2;
      }

      const integralVal = romPrev[decdigs - 1];
      const mu0 = 4 * Math.PI * 1e-7;
      return mu0 * Rs * integralVal;
    }

    // Three.js vars
    let scene, camera, renderer, controls, coilGroup;
    let txMeshes = [];
    let rxMeshes = [];

    function initThree() {
      const container = document.getElementById('coilVis');

      scene = new THREE.Scene();
      scene.background = new THREE.Color(0xf7f9fa);
      scene.up.set(0, 0, 1);

      camera = new THREE.PerspectiveCamera(
        60,
        container.clientWidth / container.clientHeight,
        0.1,
        1000
      );
      camera.up.set(0, 0, 1);
      camera.position.set(15, 15, 15);
      camera.lookAt(0, 0, 0);

      renderer = new THREE.WebGLRenderer({ antialias: true });
      renderer.setSize(container.clientWidth, container.clientHeight);
      container.appendChild(renderer.domElement);

      controls = new THREE.OrbitControls(camera, renderer.domElement);
      controls.enableDamping = true;
      controls.dampingFactor = 0.1;
      controls.screenSpacePanning = false;
      controls.minDistance = 2;
      controls.maxDistance = 200;

      const ambientLight = new THREE.AmbientLight(0xaaaaaa, 0.6);
      scene.add(ambientLight);
      const dirLight = new THREE.DirectionalLight(0xffffff, 0.6);
      dirLight.position.set(20, 20, 20);
      scene.add(dirLight);

      const axisLength = 5;
      const arrowThickness = 0.2;
      const arrowX = new THREE.ArrowHelper(
        new THREE.Vector3(1, 0, 0),
        new THREE.Vector3(0, 0, 0),
        axisLength,
        0xff0000,
        arrowThickness,
        arrowThickness * 0.5
      );
      scene.add(arrowX);

      const arrowY = new THREE.ArrowHelper(
        new THREE.Vector3(0, 1, 0),
        new THREE.Vector3(0, 0, 0),
        axisLength,
        0x00ff00,
        arrowThickness,
        arrowThickness * 0.5
      );
      scene.add(arrowY);

      const arrowZ = new THREE.ArrowHelper(
        new THREE.Vector3(0, 0, 1),
        new THREE.Vector3(0, 0, 0),
        axisLength,
        0x0000ff,
        arrowThickness,
        arrowThickness * 0.5
      );
      scene.add(arrowZ);

      addAxisLabels(axisLength);

      coilGroup = new THREE.Group();
      scene.add(coilGroup);

      window.addEventListener('resize', onWindowResize);
      animateThree();
    }

    function addAxisLabels(length) {
      function makeTextSprite(text, colorHex) {
        const canvas = document.createElement('canvas');
        canvas.width = 128;
        canvas.height = 128;
        const ctx = canvas.getContext('2d');
        ctx.fillStyle = colorHex;
        ctx.font = 'Bold 64px Arial';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.fillText(text, 64, 64);
        const texture = new THREE.CanvasTexture(canvas);
        const material = new THREE.SpriteMaterial({ map: texture, depthTest: false });
        const sprite = new THREE.Sprite(material);
        sprite.scale.set(1.5, 1.5, 1.5);
        return sprite;
      }

      const xLabel = makeTextSprite('X', '#ff0000');
      xLabel.position.set(length + 1, 0, 0);
      scene.add(xLabel);

      const yLabel = makeTextSprite('Y', '#00ff00');
      yLabel.position.set(0, length + 1, 0);
      scene.add(yLabel);

      const zLabel = makeTextSprite('Z', '#0000ff');
      zLabel.position.set(0, 0, length + 1);
      scene.add(zLabel);
    }

    function onWindowResize() {
      const container = document.getElementById('coilVis');
      camera.aspect = container.clientWidth / container.clientHeight;
      camera.updateProjectionMatrix();
      renderer.setSize(container.clientWidth, container.clientHeight);
    }

    function animateThree() {
      requestAnimationFrame(animateThree);
      controls.update();
      renderer.render(scene, camera);
    }

    function updateCoilVisuals() {
      txMeshes.forEach(m => {
        coilGroup.remove(m);
        m.geometry.dispose();
        m.material.dispose();
      });
      rxMeshes.forEach(m => {
        coilGroup.remove(m);
        m.geometry.dispose();
        m.material.dispose();
      });
      txMeshes = [];
      rxMeshes = [];

      const r_max_tx_cm = parseFloat(document.getElementById('rmaxTx').value);
      const r_max_rx_cm = parseFloat(document.getElementById('rmaxRx').value);
      const x_off_cm = parseFloat(document.getElementById('xOffset').value);
      const y_off_cm = parseFloat(document.getElementById('yOffset').value);
      const z_off_cm = parseFloat(document.getElementById('zOffset').value);
      const phiDeg = parseFloat(document.getElementById('phi').value);
      const thetaDeg = parseFloat(document.getElementById('theta').value);

      const phiRad = THREE.MathUtils.degToRad(phiDeg);
      const thetaRad = THREE.MathUtils.degToRad(thetaDeg);

      const N_spiral_tx = parseInt(document.getElementById('nSpiralTx').value);
      const N_helix_tx  = parseInt(document.getElementById('nHelixTx').value);
      const p_spiral_tx = document.getElementById('pSpiralTx').value
        .split(',').map(x => parseFloat(x.trim()) * 1e-1);
      const p_helix_tx  = document.getElementById('pHelixTx').value
        .split(',').map(x => parseFloat(x.trim()) * 1e-1);

      const N_spiral_rx = parseInt(document.getElementById('nSpiralRx').value);
      const N_helix_rx  = parseInt(document.getElementById('nHelixRx').value);
      const p_spiral_rx = document.getElementById('pSpiralRx').value
        .split(',').map(x => parseFloat(x.trim()) * 1e-1);
      const p_helix_rx  = document.getElementById('pHelixRx').value
        .split(',').map(x => parseFloat(x.trim()) * 1e-1);

      const nx = Math.sin(thetaRad) * Math.cos(phiRad);
      const ny = Math.sin(thetaRad) * Math.sin(phiRad);
      const nz = Math.cos(thetaRad);
      const nVec = new THREE.Vector3(nx, ny, nz).normalize();
      const quat = new THREE.Quaternion();
      quat.setFromUnitVectors(new THREE.Vector3(0, 0, 1), nVec);

      const tubeRadius = 0.05;

      // Transmitter: draw 모든 loop
      for (let i = 0; i < N_spiral_tx; i++) {
        const pitch_s_tx = (i < p_spiral_tx.length ? p_spiral_tx[i] : p_spiral_tx[p_spiral_tx.length - 1]);
        const radius = r_max_tx_cm - i * pitch_s_tx;
        for (let j = 0; j < N_helix_tx; j++) {
          const pitch_h_tx = (j < p_helix_tx.length ? p_helix_tx[j] : p_helix_tx[p_helix_tx.length - 1]);
          const zPos = j * pitch_h_tx;
          const geom = new THREE.TorusGeometry(radius, tubeRadius, 16, 128);
          const mat = new THREE.MeshPhongMaterial({
            color: 0x0077ff,
            shininess: 80,
            opacity: 0.7,
            transparent: true
          });
          const mesh = new THREE.Mesh(geom, mat);
          mesh.position.set(0, 0, zPos);
          coilGroup.add(mesh);
          txMeshes.push(mesh);
        }
      }

      // Receiver: draw 모든 loop
      for (let i = 0; i < N_spiral_rx; i++) {
        const pitch_s_rx = (i < p_spiral_rx.length ? p_spiral_rx[i] : p_spiral_rx[p_spiral_rx.length - 1]);
        const radius = r_max_rx_cm - i * pitch_s_rx;
        for (let j = 0; j < N_helix_rx; j++) {
          const pitch_h_rx = (j < p_helix_rx.length ? p_helix_rx[j] : p_helix_rx[p_helix_rx.length - 1]);
          const zPos = z_off_cm + j * pitch_h_rx;
          const geom = new THREE.TorusGeometry(radius, tubeRadius, 16, 128);
          const mat = new THREE.MeshPhongMaterial({
            color: 0xff3333,
            shininess: 80,
            opacity: 0.7,
            transparent: true
          });
          const mesh = new THREE.Mesh(geom, mat);
          mesh.setRotationFromQuaternion(quat);
          mesh.position.set(x_off_cm, y_off_cm, zPos);
          coilGroup.add(mesh);
          rxMeshes.push(mesh);
        }
      }
    }

    // Sweep Data 저장용
    let sweepLabels = [];
    let sweepValues = [];

    function calculateMutualInductance() {
      const r_max_tx = parseFloat(document.getElementById('rmaxTx').value) * 1e-2;
      const r_max_rx = parseFloat(document.getElementById('rmaxRx').value) * 1e-2;

      const p_spiral_tx = document.getElementById('pSpiralTx').value.split(',')
        .map(x => parseFloat(x.trim()) * 1e-3);
      const p_helix_tx = document.getElementById('pHelixTx').value.split(',')
        .map(x => parseFloat(x.trim()) * 1e-3);
      const N_spiral_tx = parseInt(document.getElementById('nSpiralTx').value);
      const N_helix_tx  = parseInt(document.getElementById('nHelixTx').value);

      const p_spiral_rx = document.getElementById('pSpiralRx').value.split(',')
        .map(x => parseFloat(x.trim()) * 1e-3);
      const p_helix_rx = document.getElementById('pHelixRx').value.split(',')
        .map(x => parseFloat(x.trim()) * 1e-3);
      const N_spiral_rx = parseInt(document.getElementById('nSpiralRx').value);
      const N_helix_rx  = parseInt(document.getElementById('nHelixRx').value);

      let x_off = parseFloat(document.getElementById('xOffset').value) * 1e-2;
      let y_off = parseFloat(document.getElementById('yOffset').value) * 1e-2;
      let z_off = parseFloat(document.getElementById('zOffset').value) * 1e-2;

      let phi   = parseFloat(document.getElementById('phi').value) * Math.PI / 180;
      let theta = parseFloat(document.getElementById('theta').value) * Math.PI / 180;

      const a = Math.sin(phi) * Math.sin(theta);
      const b = -Math.cos(phi) * Math.sin(theta);
      const c = Math.cos(theta);

      const N_tx = N_spiral_tx * N_helix_tx;
      let Rp_tx = new Array(N_tx);
      let z_tx  = new Array(N_tx);
      for (let i = 0; i < N_spiral_tx; i++) {
        const pitch_s = (i < p_spiral_tx.length ? p_spiral_tx[i] : p_spiral_tx[p_spiral_tx.length - 1]);
        const radius = r_max_tx - i * pitch_s;
        for (let j = 0; j < N_helix_tx; j++) {
          const idx = i * N_helix_tx + j;
          const pitch_h = (j < p_helix_tx.length ? p_helix_tx[j] : p_helix_tx[p_helix_tx.length - 1]);
          const zPos = j * pitch_h;
          Rp_tx[idx] = radius;
          z_tx[idx]  = zPos;
        }
      }

      const N_rx = N_spiral_rx * N_helix_rx;
      let Rs_rx = new Array(N_rx);
      let z_rx  = new Array(N_rx);
      for (let i = 0; i < N_spiral_rx; i++) {
        const pitch_s = (i < p_spiral_rx.length ? p_spiral_rx[i] : p_spiral_rx[p_spiral_rx.length - 1]);
        const radius = r_max_rx - i * pitch_s;
        for (let j = 0; j < N_helix_rx; j++) {
          const idx = i * N_helix_rx + j;
          const pitch_h = (j < p_helix_rx.length ? p_helix_rx[j] : p_helix_rx[p_helix_rx.length - 1]);
          const zPos = z_off + j * pitch_h;
          Rs_rx[idx] = radius;
          z_rx[idx]  = zPos;
        }
      }

      let Msum = 0;
      for (let i = 0; i < N_tx; i++) {
        for (let j = 0; j < N_rx; j++) {
          let zc = z_rx[j] - z_tx[i];
          let Mij = Babic_24(Rp_tx[i], Rs_rx[j], [x_off, y_off, zc], [a, b, c]);
          if (!isFinite(Mij) || isNaN(Mij)) {
            zc += 1e-8;
            Mij = Babic_24(Rp_tx[i], Rs_rx[j], [x_off, y_off, zc], [a, b, c]);
            if (!isFinite(Mij) || isNaN(Mij)) {
              Mij = 0;
            }
          }
          Msum += Mij;
        }
      }

      const M_uH = Msum * 1e6;
      document.getElementById('result').innerHTML =
        `<p><strong>Mutual Inductance M:</strong> ${M_uH.toFixed(6)} μH</p>`;

      updateCoilVisuals();
    }

    // Async version for incremental chart updates
    async function plotSweepGraph() {
      const progressDiv = document.getElementById('progress');
      progressDiv.textContent = 'Progress: 0%';

      const r_max_tx = parseFloat(document.getElementById('rmaxTx').value) * 1e-2;
      const r_max_rx = parseFloat(document.getElementById('rmaxRx').value) * 1e-2;

      const p_spiral_tx = document.getElementById('pSpiralTx').value.split(',')
        .map(x => parseFloat(x.trim()) * 1e-3);
      const p_helix_tx = document.getElementById('pHelixTx').value.split(',')
        .map(x => parseFloat(x.trim()) * 1e-3);
      const N_spiral_tx = parseInt(document.getElementById('nSpiralTx').value);
      const N_helix_tx  = parseInt(document.getElementById('nHelixTx').value);

      const p_spiral_rx = document.getElementById('pSpiralRx').value.split(',')
        .map(x => parseFloat(x.trim()) * 1e-3);
      const p_helix_rx = document.getElementById('pHelixRx').value.split(',')
        .map(x => parseFloat(x.trim()) * 1e-3);
      const N_spiral_rx = parseInt(document.getElementById('nSpiralRx').value);
      const N_helix_rx  = parseInt(document.getElementById('nHelixRx').value);

      let x_off = parseFloat(document.getElementById('xOffset').value) * 1e-2;
      let y_off = parseFloat(document.getElementById('yOffset').value) * 1e-2;
      let z_off = parseFloat(document.getElementById('zOffset').value) * 1e-2;
      let phi   = parseFloat(document.getElementById('phi').value) * Math.PI / 180;
      let theta = parseFloat(document.getElementById('theta').value) * Math.PI / 180;

      const sweepVar = document.getElementById('sweepVar').value;
      const sweepMin = parseFloat(document.getElementById('sweepMin').value);
      const sweepMax = parseFloat(document.getElementById('sweepMax').value);
      const sweepStep = parseFloat(document.getElementById('sweepStep').value);

      const points = [];
      for (let v = sweepMin; v <= sweepMax + 1e-12; v += sweepStep) {
        points.push(parseFloat(v.toFixed(6)));
      }

      // 초기화: 차트 데이터 비우기
      sweepLabels = [];
      sweepValues = [];
      window.sweepChartInstance.data.labels = [];
      window.sweepChartInstance.data.datasets[0].data = [];
      window.sweepChartInstance.update();

      for (let k = 0; k < points.length; k++) {
        // 진행률 업데이트
        const percent = Math.round(((k + 1) / points.length) * 100);
        progressDiv.textContent = `Progress: ${percent}%`;

        switch (sweepVar) {
          case 'x': x_off = points[k] * 1e-2; break;
          case 'y': y_off = points[k] * 1e-2; break;
          case 'z': z_off = points[k] * 1e-2; break;
          case 'phi': phi = points[k] * Math.PI / 180; break;
          case 'theta': theta = points[k] * Math.PI / 180; break;
        }
        const a = Math.sin(phi) * Math.sin(theta);
        const b = -Math.cos(phi) * Math.sin(theta);
        const c = Math.cos(theta);

        // Tx 배열
        const N_tx = N_spiral_tx * N_helix_tx;
        let Rp_tx = new Array(N_tx);
        let z_tx  = new Array(N_tx);
        for (let i = 0; i < N_spiral_tx; i++) {
          const pitch_s = (i < p_spiral_tx.length ? p_spiral_tx[i] : p_spiral_tx[p_spiral_tx.length - 1]);
          const radius = r_max_tx - i * pitch_s;
          for (let j = 0; j < N_helix_tx; j++) {
            const idx = i * N_helix_tx + j;
            const pitch_h = (j < p_helix_tx.length ? p_helix_tx[j] : p_helix_tx[p_helix_tx.length - 1]);
            const zPos = j * pitch_h;
            Rp_tx[idx] = radius;
            z_tx[idx]  = zPos;
          }
        }

        // Rx 배열
        const N_rx = N_spiral_rx * N_helix_rx;
        let Rs_rx = new Array(N_rx);
        let z_rx  = new Array(N_rx);
        for (let i = 0; i < N_spiral_rx; i++) {
          const pitch_s = (i < p_spiral_rx.length ? p_spiral_rx[i] : p_spiral_rx[p_spiral_rx.length - 1]);
          const radius = r_max_rx - i * pitch_s;
          for (let j = 0; j < N_helix_rx; j++) {
            const idx = i * N_helix_rx + j;
            const pitch_h = (j < p_helix_rx.length ? p_helix_rx[j] : p_helix_rx[p_helix_rx.length - 1]);
            const zPos = z_off + j * pitch_h;
            Rs_rx[idx] = radius;
            z_rx[idx]  = zPos;
          }
        }

        // Msum 계산
        let Msum = 0;
        for (let i = 0; i < N_tx; i++) {
          for (let j = 0; j < N_rx; j++) {
            let zc = z_rx[j] - z_tx[i];
            let Mij = Babic_24(Rp_tx[i], Rs_rx[j], [x_off, y_off, zc], [a, b, c]);
            if (!isFinite(Mij) || isNaN(Mij)) {
              zc += 1e-8;
              Mij = Babic_24(Rp_tx[i], Rs_rx[j], [x_off, y_off, zc], [a, b, c]);
              if (!isFinite(Mij) || isNaN(Mij)) {
                Mij = 0;
              }
            }
            Msum += Mij;
          }
        }
        const Mval = Msum * 1e6;

        // 차트 업데이트: 레이블, 데이터 추가
        sweepLabels.push(points[k].toFixed(3));
        sweepValues.push(Mval);
        window.sweepChartInstance.data.labels = sweepLabels.slice();
        window.sweepChartInstance.data.datasets[0].data = sweepValues.slice();
        window.sweepChartInstance.update();

        // 짧게 대기해서 UI가 갱신될 시간 확보
        await new Promise(resolve => setTimeout(resolve, 0));
      }

      // 완료되면 진행률 표시 제거
      progressDiv.textContent = '';
    }

    function saveSweepData() {
      if (!sweepLabels.length) {
        alert('No sweep data to save. Please plot first.');
        return;
      }
      const lines = [['Variable', 'M(μH)']];
      for (let i = 0; i < sweepLabels.length; i++) {
        lines.push([sweepLabels[i], sweepValues[i].toFixed(6)]);
      }
      const csvContent = lines.map(e => e.join(' ')).join('\n');
      const blob = new Blob([csvContent], { type: 'text/plain' });
      const url = URL.createObjectURL(blob);
      const link = document.createElement('a');
      link.href = url;
      link.download = `sweep_data.txt`;
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
      URL.revokeObjectURL(url);
    }

    window.addEventListener('load', () => {
      initThree();
      const ctx = document.getElementById('sweepChart').getContext('2d');
      window.sweepChartInstance = new Chart(ctx, {
        type: 'line',
        data: {
          labels: [],
          datasets: [{
            label: '',
            data: [],
            borderColor: '#1abc9c',
            backgroundColor: 'rgba(26, 188, 156, 0.2)',
            fill: true,
            tension: 0.2
          }]
        },
        options: {
          scales: {
            x: { display: true, title: { display: true, text: '' } },
            y: { display: true, title: { display: true, text: '' } }
          },
          plugins: { tooltip: { enabled: true } }
        }
      });
    });

    document.getElementById('calculateBtn').addEventListener('click', calculateMutualInductance);
    document.getElementById('plotBtn').addEventListener('click', () => {
      plotSweepGraph();
    });
    document.getElementById('saveDataBtn').addEventListener('click', saveSweepData);
