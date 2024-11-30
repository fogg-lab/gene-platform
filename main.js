const { app, BrowserWindow } = require('electron');
const path = require('path');

let mainWindow;

async function createWindow() {
  const isDev = (await import('electron-is-dev')).default;

  mainWindow = new BrowserWindow({
    width: 1200,
    height: 900,
    minWidth: 1200,
    minHeight: 600,
    webPreferences: {
      nodeIntegration: true,
      contextIsolation: false,
    }
  });

  const url = isDev
    ? 'http://localhost:9328'
    : `file://${path.join(__dirname, 'dist/index.html')}`;

  console.log('Loading URL:', url);

  mainWindow.loadURL(url).catch(err => {
    console.error('Failed to load URL:', err);
  });

  mainWindow.on('closed', () => {
    mainWindow = null;
  });

  if (isDev) {
    try {
      const { default: installExtension, REACT_DEVELOPER_TOOLS } = require('electron-extension-installer');
      await installExtension(REACT_DEVELOPER_TOOLS, {
        loadExtensionOptions: {
          allowFileAccess: true,
        },
      });
    } catch (e) {
      console.error('Failed to install extension:', e);
    }
  }
}

app.on("ready", async () => {
  createWindow();
});

app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') {
    app.quit();
  }
});

app.on('activate', () => {
  if (mainWindow === null) {
    createWindow();
  }
});
