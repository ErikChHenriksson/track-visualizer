//Get audio stream from browser media
navigator.mediaDevices.getDisplayMedia({ audio: true, video: true }).then(
  stream => {
    //Create audio wave from analyser data
    var audioCtx = new (window.AudioContext || window.webkitAudioContext)();
    //var masterGain = audioCtx.createGain();
    //masterGain.connect(audioCtx.destination);
    const analyser = audioCtx.createAnalyser();
    audioCtx.createMediaStreamSource(stream).connect(analyser);
    var bufferLength = analyser.frequencyBinCount;
    const waveform = new Float32Array(bufferLength);
    const freqData = new Uint8Array(bufferLength);
    // Copy the current waveform data into a Float32Array array passed into it, namely waveform.
    analyser.getFloatTimeDomainData(waveform);

    const array32 = new Float32Array(analyser.fftSize);
    analyser.getFloatTimeDomainData(array32);
    
    updateWaveform();
    drawOscilloscope();
    drawFrequencyBars();

    function updateWaveform() {
      requestAnimationFrame(updateWaveform);
      analyser.getFloatTimeDomainData(waveform);
    }
    function drawOscilloscope() {
      requestAnimationFrame(drawOscilloscope);

      const scopeCanvas = document.getElementById("oscilloscope");
      const scopeContext = scopeCanvas.getContext("2d");

      scopeCanvas.width = waveform.length;
      scopeCanvas.height = window.innerHeight * 0.2;

      scopeContext.clearRect(0, 0, scopeCanvas.width, scopeCanvas.height);
      scopeContext.beginPath();

      for (let i = 0; i < waveform.length; i++) {
        const x = i;
        const y = (0.5 + waveform[i] / 2) * scopeCanvas.height;

        if (i == 0) {
          scopeContext.moveTo(x, y);
        } else {
          scopeContext.lineTo(x, y);
        }
      }

      scopeContext.strokeStyle = "#1db954";
      scopeContext.stroke();
    }
    function drawFrequencyBars() {
      requestAnimationFrame(drawFrequencyBars);

      const freqCanvas = document.getElementById("frequency_bars");
      const freqCtx = freqCanvas.getContext("2d");
      freqCanvas.height = window.innerHeight * 0.2;
      let canvasHeight = freqCanvas.height;
      let canvasWidth = freqCanvas.width;

      analyser.getByteFrequencyData(freqData);

      freqCtx.fillStyle = "rgb(0, 0, 0)";
      freqCtx.fillRect(0, 0, canvasWidth, canvasHeight);

      var barWidth = (canvasWidth / bufferLength) * 2.5;
      var barHeight;
      var x = 0;
      for (var i = 0; i < bufferLength; i++) {
        barHeight = freqData[i];
        freqCtx.fillStyle = "rgb(50," + (barHeight + 100) + ",50)";
        freqCtx.fillRect(x, canvasHeight - barHeight / 2, barWidth, barHeight);

        x += barWidth + 1;
      }
    }
  },
  err => console.log(err)
);

anychart.onDocumentReady(function() {
  // create data
  let tag_cloud = document.getElementById("tag_cloud");
  if(tag_cloud === null) return;
  const string = tag_cloud.getAttribute("data-lyrics");
  const title = tag_cloud.getAttribute("data-title");
  let chart = anychart.tagCloud();
  chart.data(string, {
    mode: "byWord",
    maxItems: 100,
    ignoreItems: [
      "the",
      "and",
      "he",
      "or",
      "of",
      "on",
      "to",
      "in",
      "a",
      "it",
      "its",
      "is",
      "theres",
      "there"
    ]
  });
  chart.container("tag_cloud");
  chart.angles([0]);
  chart.background().fill("#000");
  chart.colorScale(
    anychart.scales
      .ordinalColor()
      .colors([
        "#327a32",
        "#328a32",
        "#329a32",
        "#32aa32",
        "#32ba32",
        "#32ca32"
      ])
  );
  chart.draw();
});