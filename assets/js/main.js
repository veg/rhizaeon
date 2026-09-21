// RhizAeon Benchmark Compendium - Interactive Scripts

document.addEventListener("DOMContentLoaded", () => {
  initLightbox();
  initCopyButtons();
  initTableFilters();
  initScalingPlot();
});

// Lightbox Modal for Figures
function initLightbox() {
  const modal = document.getElementById("lightbox-modal");
  if (!modal) return;
  const modalImg = document.getElementById("lightbox-img");
  const modalCaption = document.getElementById("lightbox-caption");
  const closeBtn = document.getElementById("lightbox-close");

  document.querySelectorAll(".figure-wrapper, .study-figure-container").forEach((el) => {
    el.addEventListener("click", () => {
      const img = el.querySelector("img");
      if (!img) return;
      modalImg.src = img.src;
      modalCaption.textContent = img.alt || "";
      modal.classList.add("active");
    });
  });

  if (closeBtn) {
    closeBtn.addEventListener("click", () => modal.classList.remove("active"));
  }

  modal.addEventListener("click", (e) => {
    if (e.target === modal || e.target === modalImg) {
      modal.classList.remove("active");
    }
  });

  document.addEventListener("keydown", (e) => {
    if (e.key === "Escape" && modal.classList.contains("active")) {
      modal.classList.remove("active");
    }
  });
}

// Copy Code Button
function initCopyButtons() {
  document.querySelectorAll(".copy-btn, .rams-copy-pill").forEach((btn) => {
    btn.addEventListener("click", () => {
      const targetId = btn.getAttribute("data-target");
      let text = "";
      if (targetId) {
        const codeEl = document.getElementById(targetId);
        if (codeEl) text = codeEl.innerText || codeEl.textContent;
      } else {
        const siblingCode = btn.parentElement ? btn.parentElement.querySelector("code") : null;
        if (siblingCode) text = siblingCode.innerText || siblingCode.textContent;
      }
      if (!text) return;
      navigator.clipboard.writeText(text).then(() => {
        const orig = btn.innerText;
        btn.innerText = "COPIED!";
        btn.style.color = "#10b981";
        setTimeout(() => {
          btn.innerText = orig;
          btn.style.color = "";
        }, 2000);
      });
    });
  });
}

// Table and Card Filtering and Search on Portal Index
function initTableFilters() {
  const searchInput = document.getElementById("study-search");
  const table = document.getElementById("benchmarks-table");
  const cards = document.querySelectorAll(".cohort-card");
  const categoryButtons = document.querySelectorAll(".filter-btn[data-category]");
  const statusButtons = document.querySelectorAll(".filter-btn[data-status]");
  const countDisplay = document.getElementById("filter-count");

  let activeCategory = "all";
  let activeStatus = "all";
  let searchQuery = "";

  function filterItems() {
    let visibleCount = 0;

    // Filter table rows
    if (table) {
      const rows = table.querySelectorAll("tbody tr");
      rows.forEach((row) => {
        const rowCategory = row.getAttribute("data-category") || "";
        const rowStatus = row.getAttribute("data-status") || "";
        const rowText = row.textContent.toLowerCase();

        const matchCategory = (activeCategory === "all" || rowCategory === activeCategory);
        const matchStatus = (activeStatus === "all" || rowStatus === activeStatus);
        const matchSearch = (!searchQuery || rowText.includes(searchQuery));

        if (matchCategory && matchStatus && matchSearch) {
          row.style.display = "";
          visibleCount++;
        } else {
          row.style.display = "none";
        }
      });
    }

    // Filter cards
    cards.forEach((card) => {
      const cardCategory = card.getAttribute("data-category") || "";
      const cardStatus = card.getAttribute("data-status") || "";
      const cardText = card.textContent.toLowerCase();

      const matchCategory = (activeCategory === "all" || cardCategory === activeCategory);
      const matchStatus = (activeStatus === "all" || cardStatus === activeStatus);
      const matchSearch = (!searchQuery || cardText.includes(searchQuery));

      if (matchCategory && matchStatus && matchSearch) {
        card.style.display = "";
      } else {
        card.style.display = "none";
      }
    });

    if (countDisplay) {
      countDisplay.textContent = `${visibleCount} of 10 cohorts shown`;
    }
  }

  if (searchInput) {
    searchInput.addEventListener("input", (e) => {
      searchQuery = e.target.value.toLowerCase().trim();
      filterItems();
    });
  }

  categoryButtons.forEach((btn) => {
    btn.addEventListener("click", () => {
      categoryButtons.forEach((b) => b.classList.remove("active"));
      btn.classList.add("active");
      activeCategory = btn.getAttribute("data-category");
      filterItems();
    });
  });

  statusButtons.forEach((btn) => {
    btn.addEventListener("click", () => {
      statusButtons.forEach((b) => b.classList.remove("active"));
      btn.classList.add("active");
      activeStatus = btn.getAttribute("data-status");
      filterItems();
    });
  });

  // Prevent link click in table from propagating to row click
  if (table) {
    table.querySelectorAll("tbody a").forEach((a) => {
      a.addEventListener("click", (e) => {
        e.stopPropagation();
      });
    });
  }
}

// Interactive SVG Scaling & Speedup Chart
function initScalingPlot() {
  const container = document.getElementById("scaling-plot-container");
  if (!container || typeof BENCHMARK_SCALING_DATA === "undefined") return;

  container.innerHTML = "";

  const width = 880;
  const height = 440;
  const margin = { top: 35, right: 45, bottom: 55, left: 75 };

  // Data: Complexity (N * L nt) vs Runtime (ms)
  const plotData = BENCHMARK_SCALING_DATA;

  // Log10 scales
  const minX = 3.5; // ~3,000 nt
  const maxX = 9.0; // ~1,000,000,000 nt
  const minY = -0.5; // ~0.3 ms
  const maxY = 6.5; // ~3,000,000 ms (~50 min)

  function scaleX(logVal) {
    return margin.left + ((logVal - minX) / (maxX - minX)) * (width - margin.left - margin.right);
  }

  function scaleY(logVal) {
    return height - margin.bottom - ((logVal - minY) / (maxY - minY)) * (height - margin.top - margin.bottom);
  }

  const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
  svg.setAttribute("viewBox", `0 0 ${width} ${height}`);
  svg.style.width = "100%";
  svg.style.height = "auto";
  svg.style.display = "block";

  // Grid lines
  const gridGroup = document.createElementNS("http://www.w3.org/2000/svg", "g");
  for (let p = 4; p <= 9; p++) {
    const x = scaleX(p);
    const line = document.createElementNS("http://www.w3.org/2000/svg", "line");
    line.setAttribute("x1", x);
    line.setAttribute("x2", x);
    line.setAttribute("y1", margin.top);
    line.setAttribute("y2", height - margin.bottom);
    line.setAttribute("stroke", "#f1f5f9");
    line.setAttribute("stroke-dasharray", "3,3");
    gridGroup.appendChild(line);

    const text = document.createElementNS("http://www.w3.org/2000/svg", "text");
    text.setAttribute("x", x);
    text.setAttribute("y", height - margin.bottom + 20);
    text.setAttribute("text-anchor", "middle");
    text.setAttribute("font-size", "11");
    text.setAttribute("fill", "#64748b");
    text.textContent = `10^${p}`;
    gridGroup.appendChild(text);
  }

  for (let p = 0; p <= 6; p += 2) {
    const y = scaleY(p);
    const line = document.createElementNS("http://www.w3.org/2000/svg", "line");
    line.setAttribute("x1", margin.left);
    line.setAttribute("x2", width - margin.right);
    line.setAttribute("y1", y);
    line.setAttribute("y2", y);
    line.setAttribute("stroke", "#f1f5f9");
    line.setAttribute("stroke-dasharray", "3,3");
    gridGroup.appendChild(line);

    const text = document.createElementNS("http://www.w3.org/2000/svg", "text");
    text.setAttribute("x", margin.left - 10);
    text.setAttribute("y", y + 4);
    text.setAttribute("text-anchor", "end");
    text.setAttribute("font-size", "11");
    text.setAttribute("fill", "#64748b");
    text.textContent = p === 0 ? "1 ms" : p === 3 ? "1 s" : "1,000 s";
    gridGroup.appendChild(text);
  }
  svg.appendChild(gridGroup);

  // Axis Titles
  const xTitle = document.createElementNS("http://www.w3.org/2000/svg", "text");
  xTitle.setAttribute("x", margin.left + (width - margin.left - margin.right) / 2);
  xTitle.setAttribute("y", height - 12);
  xTitle.setAttribute("text-anchor", "middle");
  xTitle.setAttribute("font-size", "12");
  xTitle.setAttribute("font-weight", "600");
  xTitle.setAttribute("fill", "#334155");
  xTitle.textContent = "Alignment Matrix Complexity (N taxa × L sites, nt)";
  svg.appendChild(xTitle);

  const yTitle = document.createElementNS("http://www.w3.org/2000/svg", "text");
  yTitle.setAttribute("transform", `rotate(-90)`);
  yTitle.setAttribute("x", -(margin.top + (height - margin.top - margin.bottom) / 2));
  yTitle.setAttribute("y", 20);
  yTitle.setAttribute("text-anchor", "middle");
  yTitle.setAttribute("font-size", "12");
  yTitle.setAttribute("font-weight", "600");
  yTitle.setAttribute("fill", "#334155");
  yTitle.textContent = "Wall-Clock Latency (Log Scale)";
  svg.appendChild(yTitle);

  // Tooltip
  let tooltip = document.getElementById("plot-tooltip");
  if (!tooltip) {
    tooltip = document.createElement("div");
    tooltip.id = "plot-tooltip";
    tooltip.style.position = "absolute";
    tooltip.style.padding = "10px 14px";
    tooltip.style.background = "#0f172a";
    tooltip.style.color = "#ffffff";
    tooltip.style.borderRadius = "6px";
    tooltip.style.fontSize = "0.82rem";
    tooltip.style.pointerEvents = "none";
    tooltip.style.display = "none";
    tooltip.style.zIndex = "100";
    tooltip.style.boxShadow = "0 6px 16px rgba(0, 0, 0, 0.25)";
    tooltip.style.maxWidth = "280px";
    tooltip.style.lineHeight = "1.4";
    document.body.appendChild(tooltip);
  }

  // Draw Points
  plotData.forEach((d) => {
    const cx = scaleX(d.log_complexity);
    const cyRhiz = scaleY(d.log_rhiz_ms);
    const cyHist = scaleY(d.log_hist_ms);

    // Connector line between RhizAeon and Historical
    const conn = document.createElementNS("http://www.w3.org/2000/svg", "line");
    conn.setAttribute("x1", cx);
    conn.setAttribute("x2", cx);
    conn.setAttribute("y1", cyRhiz);
    conn.setAttribute("y2", cyHist);
    conn.setAttribute("stroke", "#cbd5e1");
    conn.setAttribute("stroke-width", "1.5");
    conn.setAttribute("stroke-dasharray", "2,2");
    svg.appendChild(conn);

    // Historical Point (Gray square)
    const sq = document.createElementNS("http://www.w3.org/2000/svg", "rect");
    sq.setAttribute("x", cx - 4.5);
    sq.setAttribute("y", cyHist - 4.5);
    sq.setAttribute("width", 9);
    sq.setAttribute("height", 9);
    sq.setAttribute("fill", "#94a3b8");
    sq.setAttribute("rx", 1.5);
    svg.appendChild(sq);

    // RhizAeon Point (Emerald circle)
    const circle = document.createElementNS("http://www.w3.org/2000/svg", "circle");
    circle.setAttribute("cx", cx);
    circle.setAttribute("cy", cyRhiz);
    circle.setAttribute("r", 6);
    circle.setAttribute("fill", "#059669");
    circle.setAttribute("stroke", "#ffffff");
    circle.setAttribute("stroke-width", "2");
    circle.style.cursor = "pointer";
    circle.style.transition = "transform 0.15s ease";

    circle.addEventListener("mouseenter", (e) => {
      circle.setAttribute("r", 8);
      tooltip.style.display = "block";
      tooltip.innerHTML = `
        <div style="font-weight: 700; color: #34d399; margin-bottom: 3px;">${d.name}</div>
        <div style="font-size: 0.75rem; color: #94a3b8; margin-bottom: 6px;">${d.category}</div>
        <div><strong>Matrix:</strong> ${d.complexity_label}</div>
        <div><strong>RhizAeon:</strong> <span style="color: #34d399; font-weight: 600;">${d.rhiz_time}</span></div>
        <div><strong>Baseline:</strong> ${d.hist_time} (${d.hist_method})</div>
        <div style="margin-top: 5px; padding-top: 5px; border-top: 1px solid #334155; color: #fbbf24; font-weight: 700;">
          Speedup: ${d.speedup}
        </div>
      `;
    });

    circle.addEventListener("mousemove", (e) => {
      tooltip.style.left = `${e.pageX + 15}px`;
      tooltip.style.top = `${e.pageY - 25}px`;
    });

    circle.addEventListener("mouseleave", () => {
      circle.setAttribute("r", 6);
      tooltip.style.display = "none";
    });

    circle.addEventListener("click", () => {
      window.location.href = `studies/${d.id}/index.html`;
    });

    svg.appendChild(circle);
  });

  // Legend
  const legGroup = document.createElementNS("http://www.w3.org/2000/svg", "g");
  legGroup.setAttribute("transform", `translate(${margin.left + 15}, ${margin.top + 10})`);

  // RhizAeon Legend Item
  const rCirc = document.createElementNS("http://www.w3.org/2000/svg", "circle");
  rCirc.setAttribute("cx", 6);
  rCirc.setAttribute("cy", 6);
  rCirc.setAttribute("r", 5);
  rCirc.setAttribute("fill", "#059669");
  legGroup.appendChild(rCirc);

  const rText = document.createElementNS("http://www.w3.org/2000/svg", "text");
  rText.setAttribute("x", 18);
  rText.setAttribute("y", 10);
  rText.setAttribute("font-size", "11");
  rText.setAttribute("font-weight", "600");
  rText.setAttribute("fill", "#1e293b");
  rText.textContent = "RhizAeon Continuous Coordinate Engine";
  legGroup.appendChild(rText);

  // Baseline Legend Item
  const bSq = document.createElementNS("http://www.w3.org/2000/svg", "rect");
  bSq.setAttribute("x", 260);
  bSq.setAttribute("y", 2);
  bSq.setAttribute("width", 8);
  bSq.setAttribute("height", 8);
  bSq.setAttribute("fill", "#94a3b8");
  bSq.setAttribute("rx", 1.5);
  legGroup.appendChild(bSq);

  const bText = document.createElementNS("http://www.w3.org/2000/svg", "text");
  bText.setAttribute("x", 276);
  bText.setAttribute("y", 10);
  bText.setAttribute("font-size", "11");
  bText.setAttribute("font-weight", "500");
  bText.setAttribute("fill", "#64748b");
  bText.textContent = "Historical Baseline (RDP5 / 3SEQ / DSBM / CFML)";
  legGroup.appendChild(bText);

  svg.appendChild(legGroup);
  container.appendChild(svg);
}
