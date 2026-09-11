(function registerApi(namespace) {
    "use strict";

    namespace.api = Object.freeze({
        request(url, options = {}) {
            return window.fetch(url, {
                credentials: "same-origin",
                ...options,
            });
        },
    });
})(window.ArianeFrontend = window.ArianeFrontend || {});
