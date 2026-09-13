(function installFrontendStartupGuard() {
    "use strict";

    function revealStartupError() {
        if (document.documentElement.dataset.arianeReady === "true") return;
        const message = document.getElementById("frontend-startup-error");
        if (message) message.hidden = false;
    }

    function scheduleCheck() {
        window.setTimeout(revealStartupError, 4000);
    }

    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", scheduleCheck, { once: true });
    } else {
        scheduleCheck();
    }
})();
