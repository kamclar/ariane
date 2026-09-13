(function registerComposition(namespace) {
    "use strict";

    namespace.composeStrict = function composeStrict(contributors) {
        const result = {};
        const owners = new Map();

        for (const contributor of contributors) {
            const name = String(contributor?.name || "unnamed contributor");
            const value = contributor?.value;
            if (!value || typeof value !== "object" || Array.isArray(value)) {
                throw new TypeError(`ARIANE frontend contributor ${name} must be an object`);
            }
            for (const [key, member] of Object.entries(value)) {
                if (owners.has(key)) {
                    throw new Error(
                        `ARIANE frontend key collision: ${key} is defined by ` +
                        `${owners.get(key)} and ${name}`
                    );
                }
                owners.set(key, name);
                result[key] = member;
            }
        }
        return result;
    };
})(window.ArianeFrontend = window.ArianeFrontend || {});
