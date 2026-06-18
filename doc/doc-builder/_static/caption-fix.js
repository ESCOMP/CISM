(function () {
    function moveCaptions() {
        document.querySelectorAll("table").forEach(function (table) {
            var caption = Array.prototype.find.call(table.children, function (c) {
                return c.tagName === "CAPTION";
            });
            if (!caption) return;
            // insert above the wrapper if one exists, else above the table itself
            var target = table.closest(".wy-table-responsive") || table;
            var ext = document.createElement("div");
            ext.className = "table-caption-external";
            ext.innerHTML = caption.innerHTML;
            target.parentNode.insertBefore(ext, target);
            caption.remove();
        });
    }
    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", moveCaptions);
    } else {
        moveCaptions();
    }
})();
