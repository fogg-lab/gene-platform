function isElectron() {
    // Check if the userAgent contains 'Electron'
    if (typeof navigator !== 'undefined' && navigator.userAgent.toLowerCase().indexOf('electron') > -1) {
        return true;
    }

    // Check if the process.versions['electron'] exists
    if (typeof process !== 'undefined' && process.versions && process.versions.electron) {
        return true;
    }

    // Check if the window.electron object exists
    if (typeof window !== 'undefined' && window.electron) {
        return true;
    }

    // Not running in Electron, assume running in browser
    return false;
}

function getPublicUrl() {
    return isElectron() ? '' : process.env.PUBLIC_URL;
}

export { isElectron, getPublicUrl };
