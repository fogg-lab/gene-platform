const path = require('path');
const webpack = require('webpack');
const HtmlWebpackPlugin = require('html-webpack-plugin');
const { PyodidePlugin } = require("@pyodide/webpack-plugin");
const CopyWebpackPlugin = require('copy-webpack-plugin');
const Dotenv = require('dotenv-webpack');

module.exports = {
  entry: './src/index.js',
  output: {
    path: path.resolve(__dirname, 'dist'),
    filename: 'bundle.js',
    publicPath: '/',
    globalObject: 'self'
  },
  mode: 'development',
  target: 'web',
  module: {
    rules: [
      {
        test: /\.(js|jsx)$/,
        exclude: /node_modules/,
        use: 'babel-loader'
      },
      {
        test: /\.css$/,
        use: ['style-loader', 'css-loader', 'postcss-loader'],
      },
      {
        test: /\.(woff|woff2|eot|ttf|otf)$/i,
        type: 'asset/resource',
        generator: {
          filename: 'assets/fonts/[name][ext]'
        }
      },
      {
        test: /\.(png|svg|jpg|jpeg|gif)$/i,
        use: [{
          loader: 'file-loader',
          options: {
            name: '[name].[ext]',
            outputPath: 'assets/',
            publicPath: 'assets/'
          }
        }]
      },
      {
        test: /\.worker\.js$/,
        use: {
          loader: 'worker-loader',
          options: {
            filename: '[name].[contenthash].worker.js',
            publicPath: '/'
          }
        }
      }
    ]
  },
  resolve: {
    extensions: [".js", ".jsx"],
    alias: {
      global: path.resolve(__dirname, "./global-shim.js"),
    },
    fallback: {
      "buffer": require.resolve("buffer/"),
      "stream": require.resolve("stream-browserify"),
      "process": require.resolve("process/browser"),
      "path": require.resolve("path-browserify"),
      "crypto": require.resolve("crypto-browserify"),
      "events": require.resolve("events/"),
      "util": require.resolve("util/"),
      "assert": require.resolve("assert/"),
      "fs": false,
      "net": false,
      "tls": false
    }
  },
  plugins: [
    new webpack.ProvidePlugin({
      Buffer: ['buffer', 'Buffer'],
      process: 'process/browser'
    }),
    new HtmlWebpackPlugin({
      template: './public/index.html',
      filename: 'index.html'
    }),
    new webpack.HotModuleReplacementPlugin(),
    new PyodidePlugin(),
    new CopyWebpackPlugin({
      patterns: [
        { from: 'public/py-wheels', to: 'py-wheels' }
      ],
    }),
    new Dotenv(),
    new webpack.DefinePlugin({
      'process.type': '"renderer"'
    })
  ],
  devServer: {
    static: {
      directory: path.join(__dirname, 'public')
    },
    compress: true,
    port: 9328,
    historyApiFallback: true,
    proxy: {
      '/datasets': {
        target: 'https://docgl1or94tw4.cloudfront.net',
        changeOrigin: true,
        pathRewrite: { '^/datasets': '' },
      },
      '/broadinstitute': {
        target: 'https://data.broadinstitute.org',
        changeOrigin: true,
        pathRewrite: { '^/broadinstitute': '' },
      },
      '/msigdb': {
        target: 'https://www.gsea-msigdb.org',
        changeOrigin: true,
        pathRewrite: { '^/msigdb': '' },
      },
    },
  }
};
