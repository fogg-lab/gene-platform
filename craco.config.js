const { PyodidePlugin } = require('@pyodide/webpack-plugin');
const CopyWebpackPlugin = require('copy-webpack-plugin');
const webpack = require('webpack');
const path = require('path');

module.exports = {
  webpack: {
    configure: (webpackConfig, { env, paths }) => {
      // Add proper worker-loader configuration
      webpackConfig.module.rules.push({
        test: /\.worker\.(js|ts)$/,
        use: {
          loader: 'worker-loader',
          options: {
            inline: 'no-fallback'
          }
        }
      });

      // Add fallbacks for web environment
      webpackConfig.resolve.fallback = {
        ...webpackConfig.resolve.fallback,
        "buffer": require.resolve("buffer"),
        "stream": require.resolve("stream-browserify"),
        "process": require.resolve("process/browser"),
        "assert": require.resolve("assert/"),
        "crypto": require.resolve("crypto-browserify"),
        "fs": false,
        "path": false,
        "electron": false
      };

      // Add ProvidePlugin
      webpackConfig.plugins.push(
        new webpack.ProvidePlugin({
          Buffer: ['buffer', 'Buffer'],
          process: 'process/browser',
        })
      );

      // Add PyodidePlugin
      webpackConfig.plugins.push(new PyodidePlugin());

      // Add CopyWebpackPlugin with modified configuration
      webpackConfig.plugins.push(
        new CopyWebpackPlugin({
          patterns: [
            { 
              from: path.resolve(__dirname, 'public/py-wheels'),
              to: 'static/py-wheels'
            },
            {
              from: path.resolve(__dirname, 'public/404.html'),
              to: '404.html'
            }
          ],
        })
      );

      // Define PUBLIC_URL
      const publicPath = process.env.NODE_ENV === 'production' 
        ? '/gene-platform/'
        : '/';
        
      webpackConfig.plugins.push(
        new webpack.DefinePlugin({
          'process.env.PUBLIC_URL': JSON.stringify(publicPath)
        })
      );

      // Update output configuration
      webpackConfig.output = {
        ...webpackConfig.output,
        globalObject: 'self',
        publicPath,
        filename: 'static/js/[name].[contenthash:8].js',
        chunkFilename: 'static/js/[name].[contenthash:8].chunk.js'
      };

      return webpackConfig;
    },
  },
  devServer: {
    headers: {
      "Access-Control-Allow-Origin": "*",
      "Access-Control-Allow-Methods": "GET, POST, PUT, DELETE, PATCH, OPTIONS",
      "Access-Control-Allow-Headers": "X-Requested-With, content-type, Authorization"
    },
    proxy: {
      '/gene-platform/datasets': {
        target: 'https://docgl1or94tw4.cloudfront.net',
        changeOrigin: true,
        pathRewrite: { '^/gene-platform/datasets': '' },
      },
      '/gene-platform/broadinstitute': {
        target: 'https://data.broadinstitute.org',
        changeOrigin: true,
        pathRewrite: { '^/gene-platform/broadinstitute': '' },
      },
      '/gene-platform/msigdb': {
        target: 'https://www.gsea-msigdb.org',
        changeOrigin: true,
        pathRewrite: { '^/gene-platform/msigdb': '' },
      },
    },
  }
};
