require('dotenv').config();

const { app, BrowserWindow, ipcMain } = require('electron');
const path = require('path');

let mainWindow;

ipcMain.handle('getEnvVar', (event, varName) => {
    return process.env[varName];
});

async function createWindow() {
  const isDev = (await import('electron-is-dev')).default;

  mainWindow = new BrowserWindow({
    width: 1200,
    height: 900,
    minWidth: 1200,
    minHeight: 600,
    webPreferences: {
      nodeIntegration: false,
      contextIsolation: true,
      preload: path.join(__dirname, 'preload.js')
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
